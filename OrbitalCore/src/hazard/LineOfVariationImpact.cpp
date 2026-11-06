#include "orbital_core/hazard/LineOfVariationImpact.h"
#include <functional>
#include <utility>

#include "NBodyPropagator.h"
#include "hazard/LineOfVariation.h"
#include "hazard/LovCoordinateChart.h"
#include "hazard/HazardAudit.h"
#include "hazard/lov/LovTraceDetail.h"
#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/hazard/PlanetData.h"
#include "logger/Logger.h"

#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <iomanip>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace orbital_core::hazard {

	using namespace lov_detail;

	namespace {
#ifdef _WIN32
		using WorkerAffinityMask = DWORD_PTR;

		[[nodiscard]] std::vector<WorkerAffinityMask> workerAffinityMasks()
		{
			DWORD_PTR processMask = 0;
			DWORD_PTR systemMask = 0;
			if (!GetProcessAffinityMask(GetCurrentProcess(), &processMask, &systemMask) || processMask == 0) {
				return {};
			}

			DWORD bytes = 0;
			GetLogicalProcessorInformation(nullptr, &bytes);
			if (GetLastError() != ERROR_INSUFFICIENT_BUFFER || bytes == 0)
				return {};

			std::vector<unsigned char> storage(bytes);
			auto* topology = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION>(storage.data());
			if (!GetLogicalProcessorInformation(topology, &bytes))
				return {};

			const std::size_t count = bytes / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
			std::vector<WorkerAffinityMask> primary;
			std::vector<WorkerAffinityMask> siblings;
			for (std::size_t i = 0; i < count; ++i) {
				if (topology[i].Relationship != RelationProcessorCore)
					continue;

				WorkerAffinityMask mask = topology[i].ProcessorMask & processMask;
				if (mask == 0)
					continue;

				const WorkerAffinityMask first = mask & (~mask + 1);
				primary.push_back(first);
				mask &= ~first;
				while (mask != 0) {
					const WorkerAffinityMask bit = mask & (~mask + 1);
					siblings.push_back(bit);
					mask &= ~bit;
				}
			}

			primary.insert(primary.end(), siblings.begin(), siblings.end());
			if (!primary.empty())
				return primary;

			while (processMask != 0) {
				const WorkerAffinityMask bit = processMask & (~processMask + 1);
				primary.push_back(bit);
				processMask &= ~bit;
			}
			return primary;
		}

		class ScopedWorkerAffinity final {
		   public:
			ScopedWorkerAffinity(const std::vector<WorkerAffinityMask>& masks, int workerIndex) noexcept
			{
				if (workerIndex < 0 || static_cast<std::size_t>(workerIndex) >= masks.size()) {
					return;
				}
				m_previous = SetThreadAffinityMask(GetCurrentThread(), masks[static_cast<std::size_t>(workerIndex)]);
			}

			~ScopedWorkerAffinity()
			{
				if (m_previous != 0) {
					SetThreadAffinityMask(GetCurrentThread(), m_previous);
				}
			}

			ScopedWorkerAffinity(const ScopedWorkerAffinity&) = delete;
			ScopedWorkerAffinity& operator=(const ScopedWorkerAffinity&) = delete;

		   private:
			WorkerAffinityMask m_previous{0};
		};
#else
		using WorkerAffinityMask = std::uintptr_t;

		[[nodiscard]] std::vector<WorkerAffinityMask> workerAffinityMasks() { return {}; }

		class ScopedWorkerAffinity final {
		   public:
			ScopedWorkerAffinity(const std::vector<WorkerAffinityMask>&, int) noexcept {}
		};
#endif

		struct TimingMetricSnapshot {
			std::uint64_t calls{};
			std::uint64_t nanoseconds{};
		};

		struct TimingSnapshot {
			TimingMetricSnapshot closeApproachScan;
			TimingMetricSnapshot targetPlaneDifferential;
			TimingMetricSnapshot minimumSearch;
			TimingMetricSnapshot offLovRefinement;
			std::uint64_t stmSegments{};
			std::uint64_t stmAcceptedSteps{};
			std::uint64_t stmRejectedSteps{};
			std::uint64_t stmRhsEvaluations{};
			std::uint64_t stmCorrectorIterations{};
			std::uint64_t stmIntegratedMicrodays{};
		};

		[[nodiscard]] TimingMetricSnapshot snapshotMetric(const LovTimingMetric& metric) noexcept
		{
			return {metric.calls.load(std::memory_order_relaxed), metric.nanoseconds.load(std::memory_order_relaxed)};
		}

		[[nodiscard]] TimingSnapshot snapshotTimings(const LovTimingCounters* timings) noexcept
		{
			if (!timings)
				return {};
			TimingSnapshot snapshot;
			snapshot.closeApproachScan = snapshotMetric(timings->closeApproachScan);
			snapshot.targetPlaneDifferential = snapshotMetric(timings->targetPlaneDifferential);
			snapshot.minimumSearch = snapshotMetric(timings->minimumSearch);
			snapshot.offLovRefinement = snapshotMetric(timings->offLovRefinement);
			snapshot.stmSegments = timings->stmSegments.load(std::memory_order_relaxed);
			snapshot.stmAcceptedSteps = timings->stmAcceptedSteps.load(std::memory_order_relaxed);
			snapshot.stmRejectedSteps = timings->stmRejectedSteps.load(std::memory_order_relaxed);
			snapshot.stmRhsEvaluations = timings->stmRhsEvaluations.load(std::memory_order_relaxed);
			snapshot.stmCorrectorIterations = timings->stmCorrectorIterations.load(std::memory_order_relaxed);
			snapshot.stmIntegratedMicrodays = timings->stmIntegratedMicrodays.load(std::memory_order_relaxed);
			return snapshot;
		}

		[[nodiscard]] std::uint64_t deltaCalls(const TimingMetricSnapshot& current,
		                                       const TimingMetricSnapshot& baseline) noexcept
		{
			return current.calls - baseline.calls;
		}

		[[nodiscard]] double deltaSeconds(const TimingMetricSnapshot& current,
		                                  const TimingMetricSnapshot& baseline) noexcept
		{
			return static_cast<double>(current.nanoseconds - baseline.nanoseconds) * 1.0e-9;
		}

		class ScopedMetricTimer final {
		   public:
			explicit ScopedMetricTimer(LovTimingMetric* metric) noexcept
			    : m_metric(metric), m_started(std::chrono::steady_clock::now())
			{
			}

			~ScopedMetricTimer()
			{
				if (!m_metric)
					return;
				const auto elapsed =
				    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - m_started)
				        .count();
				m_metric->nanoseconds.fetch_add(static_cast<std::uint64_t>(elapsed), std::memory_order_relaxed);
				m_metric->calls.fetch_add(1, std::memory_order_relaxed);
			}

		   private:
			LovTimingMetric* m_metric;
			std::chrono::steady_clock::time_point m_started;
		};

		class PhaseProgressLogger final {
		   public:
			PhaseProgressLogger(std::string phase, int total, const std::atomic<int>& done,
			                    const std::atomic<int>* active = nullptr,
			                    const std::atomic<std::uint64_t>* workItems = nullptr,
			                    const LovTimingCounters* timings = nullptr)
			    : m_phase(std::move(phase)),
			      m_total(total),
			      m_done(done),
			      m_active(active),
			      m_workItems(workItems),
			      m_timings(timings),
			      m_timingBaseline(snapshotTimings(timings)),
			      m_ephemerisBaseline(DE440Source::getCacheStats()),
			      m_started(std::chrono::steady_clock::now()),
			      m_thread([this] { run(); })
			{
				logSnapshot("start");
			}

			~PhaseProgressLogger() { stop(); }

			PhaseProgressLogger(const PhaseProgressLogger&) = delete;
			PhaseProgressLogger& operator=(const PhaseProgressLogger&) = delete;

			void finish()
			{
				stop();
				logSnapshot("done");
			}

		   private:
			void run()
			{
				std::unique_lock lock(m_mutex);
				while (!m_stop) {
					if (m_wakeup.wait_for(lock, std::chrono::seconds(10), [this] { return m_stop; })) {
						break;
					}
					lock.unlock();
					logSnapshot("running");
					lock.lock();
				}
			}

			void stop()
			{
				if (!m_thread.joinable())
					return;
				{
					std::lock_guard lock(m_mutex);
					m_stop = true;
				}
				m_wakeup.notify_one();
				m_thread.join();
			}

			void logSnapshot(const char* state) const
			{
				const int done = m_done.load(std::memory_order_relaxed);
				const double elapsed =
				    std::chrono::duration<double>(std::chrono::steady_clock::now() - m_started).count();
				const double rate = elapsed > 0.0 ? done / elapsed : 0.0;
				const double eta = rate > 0.0 ? std::max(0, m_total - done) / rate : 0.0;

				auto line = LOG_INFO;
				line << std::fixed << std::setprecision(1) << "[LineOfVariation] PROGRESS"
				     << "  phase=" << m_phase << "  state=" << state << "  done=" << done << '/' << m_total
				     << "  percent=" << (m_total > 0 ? 100.0 * done / m_total : 100.0) << "  ratePerSec=" << rate;
				if (done > 0 && done < m_total)
					line << "  etaSec=" << eta;
				if (m_active)
					line << "  active=" << m_active->load(std::memory_order_relaxed);
				if (m_workItems)
					line << "  records=" << m_workItems->load(std::memory_order_relaxed);
				if (m_timings) {
					const TimingSnapshot current = snapshotTimings(m_timings);
					const std::uint64_t stmSegments = current.stmSegments - m_timingBaseline.stmSegments;
					const std::uint64_t accepted = current.stmAcceptedSteps - m_timingBaseline.stmAcceptedSteps;
					const std::uint64_t rejected = current.stmRejectedSteps - m_timingBaseline.stmRejectedSteps;
					const std::uint64_t rhs = current.stmRhsEvaluations - m_timingBaseline.stmRhsEvaluations;
					const std::uint64_t corrector =
					    current.stmCorrectorIterations - m_timingBaseline.stmCorrectorIterations;
					const double integratedDays = 1.0e-6 * static_cast<double>(current.stmIntegratedMicrodays -
					                                                           m_timingBaseline.stmIntegratedMicrodays);
					line << "  caScanCalls="
					     << deltaCalls(current.closeApproachScan, m_timingBaseline.closeApproachScan)
					     << "  caScanWorkerSec="
					     << deltaSeconds(current.closeApproachScan, m_timingBaseline.closeApproachScan)
					     << "  tpStmCalls="
					     << deltaCalls(current.targetPlaneDifferential, m_timingBaseline.targetPlaneDifferential)
					     << "  tpStmWorkerSec="
					     << deltaSeconds(current.targetPlaneDifferential, m_timingBaseline.targetPlaneDifferential)
					     << "  minimumWorkerSec=" << deltaSeconds(current.minimumSearch, m_timingBaseline.minimumSearch)
					     << "  offLovWorkerSec="
					     << deltaSeconds(current.offLovRefinement, m_timingBaseline.offLovRefinement)
					     << "  stmSegments=" << stmSegments << "  stmAccepted=" << accepted
					     << "  stmRejected=" << rejected << "  stmRhs=" << rhs << "  stmCorrectorIters=" << corrector;
					if (accepted > 0) {
						line << "  meanStepDays=" << integratedDays / accepted
						     << "  rhsPerStep=" << static_cast<double>(rhs) / accepted
						     << "  correctorPerStep=" << static_cast<double>(corrector) / accepted;
					}

					const DE440Source::CacheStats eph = DE440Source::getCacheStats();
					const std::uint64_t ephHits = eph.hits - m_ephemerisBaseline.hits;
					const std::uint64_t ephMisses = eph.misses - m_ephemerisBaseline.misses;
					const std::uint64_t ephLookups = ephHits + ephMisses;
					line << "  ephHits=" << ephHits << "  ephMisses=" << ephMisses
					     << "  ephHitPct=" << (ephLookups > 0 ? 100.0 * ephHits / ephLookups : 0.0)
					     << "  jplStateWorkerSec="
					     << 1.0e-9 *
					            static_cast<double>(eph.jplStateNanoseconds - m_ephemerisBaseline.jplStateNanoseconds);
				}
			}

			std::string m_phase;
			int m_total;
			const std::atomic<int>& m_done;
			const std::atomic<int>* m_active;
			const std::atomic<std::uint64_t>* m_workItems;
			const LovTimingCounters* m_timings;
			TimingSnapshot m_timingBaseline;
			DE440Source::CacheStats m_ephemerisBaseline;
			std::chrono::steady_clock::time_point m_started;
			std::mutex m_mutex;
			std::condition_variable m_wakeup;
			bool m_stop{false};
			std::thread m_thread;
		};
	}  // namespace

	std::vector<ImpactRiskRow> traceLineOfVariationImpacts(const Orbit& nominalOrbit, const StateCov6x6& cov0,
	                                                       const NBodyPropagator& propagator,
	                                                       const LineOfVariationImpactConfig& cfg,
	                                                       IEphemerisSource* eph,
	                                                       std::function<bool(int, int)> progressCb)
	{
		FitParameterCovariance cov(6, 6);
		cov = cov0;
		return traceLineOfVariationImpacts(nominalOrbit, cov, propagator.nonGrav(), propagator, cfg, eph,
		                                   std::move(progressCb));
	}

	std::vector<ImpactRiskRow> traceLineOfVariationImpacts(
	    const Orbit& nominalOrbit, const FitParameterCovariance& parameterCovariance,
	    const NonGravModel& nominalNonGrav, const NBodyPropagator& propagator, const LineOfVariationImpactConfig& cfg,
	    IEphemerisSource* eph, std::function<bool(int, int)> progressCb)
	{
		std::vector<ImpactRiskRow> rows;
		const VecX xNom = nominalParameters(nominalOrbit, nominalNonGrav);
		if (!validConfig(cfg) || parameterCovariance.rows() != xNom.size() ||
		    parameterCovariance.cols() != xNom.size() || parameterCovariance.isZero(0.0)) {
			return rows;
		}

		if (!eph)
			eph = getDefaultEphemerisSource();
		if (!eph)
			return rows;

		LovCoordinateChart chart =
		    makeAutomaticLovChart(nominalOrbit, parameterCovariance, nominalNonGrav, cfg.coordinateSystem);

		if (!chart.valid) {
			LOG_INFO << "[LineOfVariation] skipped: coordinate chart unavailable";
			return rows;
		}

		const VecX units = lovCoordinateUnits(chart.coordinateSystem, chart.nominalCoordinates, cfg.coordinateScaling);
		if (units.size() != chart.nominalCoordinates.size())
			return rows;
		const VecX scales = units.cwiseInverse();
		const LovWeakDirection weak =
		    weakDirectionFromCovariance(scales.asDiagonal() * chart.coordinateCovariance * scales.asDiagonal());
		if (!weak.valid || weak.parameterStepPerSigma.size() != chart.nominalCoordinates.size()) {
			return rows;
		}

		LOG_INFO << "[LineOfVariation] coordinates=" << lovCoordinateCode(chart.coordinateSystem)
		         << " scalingLov=" << (cfg.coordinateScaling ? "on" : "off");

		LovTimingCounters timingCounters;
		const NBodyPropagator encounterPropagator =
		    propagator.forEncounterMonitoring(*eph, cfg.encounterStepPoints, cfg.closeApproach.encounterRadiusAu);

		TraceEvaluationContext traceCtx{
		    chart, parameterCovariance, nominalNonGrav, nominalOrbit, encounterPropagator, cfg, *eph};
		traceCtx.timings = &timingCounters;

		LovSamplingConfig sampling;
		sampling.sigmaMax = cfg.sigmaMax;
		sampling.sigmaStepMax = cfg.sigmaStepMax;
		sampling.probabilityStep = cfg.probabilityStep;
		sampling.probabilitySampling = cfg.probabilitySampling;
		sampling.diskRadiusBodyRadii = cfg.closeApproach.thresholdAu * kAuKm / planetRadiusKm(cfg.closeApproach.planet);
		sampling.maxPositiveSamples = cfg.maxPositiveSamples;

		const std::vector<double> sigmas = sampleLovSigmas(sampling);
		if (sigmas.empty())
			return rows;

#ifdef _OPENMP
		omp_set_dynamic(0);
		const int lovWorkerCount = defaultLovWorkerCount(omp_get_max_threads());
#else
		constexpr int lovWorkerCount = 1;
#endif
		const std::vector<WorkerAffinityMask> affinityMasks = workerAffinityMasks();

		LOG_INFO << "[LineOfVariation] TIMING"
		         << "  phase=VA-catalog"
		         << "  state=start"
		         << "  samples=" << sigmas.size()
		         << "  sampling=" << (cfg.probabilitySampling ? "probability" : "fixed")
		         << "  sigmaStepMax=" << cfg.sigmaStepMax << "  dmeaAu=" << cfg.closeApproach.encounterRadiusAu
		         << "  dterAu=" << Constants::EARTH_LOCAL_DYNAMICS_AU
		         << "  npoint=" << encounterPropagator.encounterStepPointCount();
		const auto catalogStarted = std::chrono::steady_clock::now();
		std::vector<LovCatalogPoint> catalog = buildLovCatalog(traceCtx, sigmas);
		const double catalogElapsed =
		    std::chrono::duration<double>(std::chrono::steady_clock::now() - catalogStarted).count();
		const std::size_t catalogValid = static_cast<std::size_t>(
		    std::count_if(catalog.begin(), catalog.end(), [](const LovCatalogPoint& point) { return point.valid; }));
		LOG_INFO << std::fixed << std::setprecision(3) << "[LineOfVariation] TIMING"
		         << "  phase=VA-catalog"
		         << "  state=done"
		         << "  wallSec=" << catalogElapsed << "  valid=" << catalogValid << '/' << catalog.size();

		const char* catalogMode = cfg.impactRefinement.rmsScaling.linearization ? "rms-constrained" : "linear";
		traceCtx.catalog = &catalog;

		std::vector<int> validCatalogIndices;
		validCatalogIndices.reserve(catalogValid);
		for (std::size_t i = 0; i < catalog.size(); ++i) {
			if (catalog[i].valid)
				validCatalogIndices.push_back(static_cast<int>(i));
		}

		LOG_INFO << "[LineOfVariation] scheduler"
		         << "  workers=" << lovWorkerCount
		         << "  affinity=" << (affinityMasks.empty() ? "operating-system" : "pinned-physical-first")
		         << "  affinitySlots=" << affinityMasks.size();

		std::vector<std::vector<TraceRecord>> buckets(sigmas.size());
		std::atomic<int> vaDone{0};
		std::atomic<std::uint64_t> coarseRecords{0};
		PhaseProgressLogger vaProgress("VA-propagation", static_cast<int>(validCatalogIndices.size()), vaDone, nullptr,
		                               &coarseRecords, &timingCounters);

#pragma omp parallel num_threads(lovWorkerCount)
		{
#ifdef _OPENMP
			ScopedWorkerAffinity workerAffinity(affinityMasks, omp_get_thread_num());
#else
			ScopedWorkerAffinity workerAffinity(affinityMasks, 0);
#endif
#pragma omp for schedule(dynamic, 1)
			for (int workIndex = 0; workIndex < static_cast<int>(validCatalogIndices.size()); ++workIndex) {
				const int i = validCatalogIndices[static_cast<std::size_t>(workIndex)];
				auto& bucket = buckets[static_cast<std::size_t>(i)];
				bucket = evaluateSigma(traceCtx, catalog[static_cast<std::size_t>(i)]);
				coarseRecords.fetch_add(bucket.size(), std::memory_order_relaxed);
				vaDone.fetch_add(1, std::memory_order_relaxed);
			}
		}
		vaProgress.finish();

		std::vector<TraceRecord> records;
		for (auto& bucket : buckets) {
			for (auto& record : bucket)
				records.push_back(std::move(record));
		}

		std::vector<std::vector<TraceRecord>> returns = buildReturnRecords(std::move(records));

		LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] trace"
		         << "  sigmaSamples=" << sigmas.size() << "  returns=" << returns.size() << "  maxRows=" << cfg.maxRows
		         << "  workers=" << lovWorkerCount << "  chart=" << lovCoordinateName(chart.coordinateSystem)
		         << "  catalog=" << catalogMode << "  catalogValid=" << catalogValid << "  sigmaMax=" << cfg.sigmaMax
		         << "  sigmaStepMax=" << cfg.sigmaStepMax
		         << "  sampling=" << (cfg.probabilitySampling ? "probability" : "fixed")
		         << "  probabilityStep=" << cfg.probabilityStep << "  dmeaAu=" << cfg.closeApproach.encounterRadiusAu
		         << "  npoint=" << encounterPropagator.encounterStepPointCount()
		         << "  dmeacontrAu=" << cfg.closeApproach.thresholdAu << "  effectiveRecordRadiusAu="
		         << std::min(cfg.closeApproach.encounterRadiusAu, cfg.closeApproach.thresholdAu);

		std::vector<std::vector<ImpactRiskRow>> slots(returns.size());
		std::atomic<int> progressDone{0};
		std::atomic<int> activeReturns{0};
		std::atomic<std::uint64_t> refinedRows{0};
		std::atomic<bool> cancelled{false};
		const int totalReturns = static_cast<int>(returns.size());
		PhaseProgressLogger returnProgress("return-analysis", totalReturns, progressDone, &activeReturns, &refinedRows,
		                                   &timingCounters);

		const auto markReturnDone = [&] {
			const int done = progressDone.fetch_add(1, std::memory_order_relaxed) + 1;
			if (progressCb) {
#pragma omp critical(lov_progress)
				{
					if (!progressCb(done, totalReturns))
						cancelled.store(true, std::memory_order_relaxed);
				}
			}
		};

#pragma omp parallel num_threads(lovWorkerCount)
		{
#ifdef _OPENMP
			ScopedWorkerAffinity workerAffinity(affinityMasks, omp_get_thread_num());
#else
			ScopedWorkerAffinity workerAffinity(affinityMasks, 0);
#endif
#pragma omp for schedule(dynamic, 1)
			for (int i = 0; i < totalReturns; ++i) {
				if (cancelled.load(std::memory_order_relaxed))
					continue;

				std::vector<TraceRecord>& ret = returns[static_cast<std::size_t>(i)];
				if (ret.empty()) {
					markReturnDone();
					continue;
				}

				activeReturns.fetch_add(1, std::memory_order_relaxed);

				try {
					const std::size_t rawSamples = ret.size();
					const int rawImpacts = countImpacts(ret);
					constexpr long long resampled = 0;

					std::vector<TraceRecord> denseReps = addDenseImpactSegments(ret, sigmas);
					const bool denseImpact = !denseReps.empty();

					std::vector<TraceRecord> keyholes;
					if (denseImpact) {
						keyholes = std::move(denseReps);
					} else {
						ScopedMetricTimer timing(&timingCounters.minimumSearch);
						keyholes = addRefinedMinima(ret, traceCtx);
					}

					if (keyholes.empty()) {
						LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] row"
						         << "  return=" << i << "  samples=" << ret.size() << "  rawSamples=" << rawSamples
						         << "  resampled=" << resampled << "  rawImpacts=" << rawImpacts << "  candidate=no"
						         << "  skipped=no-minimum";
						activeReturns.fetch_sub(1, std::memory_order_relaxed);
						markReturnDone();
						continue;
					}

					std::sort(keyholes.begin(), keyholes.end(),
					          [](const TraceRecord& a, const TraceRecord& b) { return a.trace.sigma < b.trace.sigma; });

					std::vector<ImpactRiskRow> returnRows;
					for (TraceRecord& rep : keyholes) {
						bool viRefined = false;
						const char* offLovCheck = "not-run";
						audit::Scope riskAudit(rep.sampleIndex, rep.trace.sigma, true);
						audit::write("risk-input", [&](auto& log) {
							log << " return=" << i << " epochJd=" << nominalOrbit.getEpoch()
							    << " denseImpact=" << (denseImpact ? "yes" : "no") << " jd=" << rep.event.jd
							    << " qAu=" << rep.event.distAu << " b=" << rep.trace.distanceBodyRadii
							    << " xiKm=" << rep.trace.coordsKm.x() << " zetaKm=" << rep.trace.coordsKm.y()
							    << " bodyRadiusKm=" << rep.trace.bodyRadiusKm
							    << " impactRadius=" << rep.trace.impactCrossSectionBodyRadii
							    << " derivative=" << rep.trace.dDistance2DSigmaBodyRadii2
							    << " stretch=" << rep.trace.stretchBodyRadii
							    << " stretchLov=" << rep.trace.stretchLovBodyRadiiPerSigma
							    << " width=" << rep.trace.widthBodyRadii << " alpha=" << rep.trace.alphaRad
							    << " gamma=" << rep.gamma << " sigmaLov=" << rep.trace.sigma
							    << " sigmaImpact=" << rep.probability.sigmaImpact
							    << " p=" << rep.probability.probability;
							for (int p = 0; p < rep.parameters.size(); ++p)
								log << " x" << p << '=' << rep.parameters(p);
							for (int p = 0; p < rep.stepPerSigma.size(); ++p)
								log << " step" << p << '=' << rep.stepPerSigma(p);
							for (int p = 0; p < rep.covariance.rows(); ++p)
								for (int q = 0; q < rep.covariance.cols(); ++q)
									log << " cov" << p << '_' << q << '=' << rep.covariance(p, q);
							for (int p = 0; p < rep.dCoordsDParams.rows(); ++p)
								for (int q = 0; q < rep.dCoordsDParams.cols(); ++q)
									log << " j" << p << '_' << q << '=' << rep.dCoordsDParams(p, q);
							for (std::size_t p = 0; p < rep.stmCheckpointJds.size(); ++p)
								log << " checkpoint" << p << '=' << rep.stmCheckpointJds[p];
						});
						if (!denseImpact) {

							if (!isReportableRisk(rep)) {
								LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] row"
								         << "  return=" << i << "  sample=" << rep.sampleIndex
								         << "  samples=" << ret.size() << "  rawSamples=" << rawSamples
								         << "  resampled=" << resampled << "  rawImpacts=" << rawImpacts
								         << "  denseImpact=no"
								         << "  viRefined=no"
								         << "  jd=" << rep.event.jd << "  missKm=" << rep.trace.distanceKm
								         << "  distanceBodyR=" << rep.trace.distanceBodyRadii
								         << "  localMoidBodyR=" << localMoidBodyRadii(rep) << "  minPossibleBodyR="
								         << minimumPossibleBodyRadii(rep, traceCtx.cfg.sigmaMax)
								         << "  widthBodyR=" << rep.trace.widthBodyRadii
								         << "  stretchBodyR=" << rep.trace.stretchBodyRadii
								         << "  stretchLovBodyR=" << rep.trace.stretchLovBodyRadiiPerSigma
								         << "  sigmaLov=" << rep.trace.sigma
								         << "  sigmaImpact=" << rep.probability.sigmaImpact << "  p=" << std::scientific
								         << std::setprecision(12) << rep.probability.probability << std::fixed
								         << std::setprecision(6)
								         << "  alphaDeg=" << rep.trace.alphaRad * Constants::RAD_TO_DEG
								         << "  alphaLovDeg=" << rep.trace.alphaLovRad * Constants::RAD_TO_DEG
								         << "  dD2dSigma=" << rep.trace.dDistance2DSigmaBodyRadii2 << "  candidate=no"
								         << "  skipped=no-risk";
								continue;
							}

							OffLovCheckResult checked;
							{
								ScopedMetricTimer timing(&timingCounters.offLovRefinement);
								checked = refineVirtualImpactorRecord(rep, traceCtx);
							}
							offLovCheck = checked.reason;
							if (checked.status == OffLovCheckStatus::Rejected) {
								LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] row"
								         << "  return=" << i << "  sample=" << rep.sampleIndex
								         << "  samples=" << ret.size() << "  rawSamples=" << rawSamples
								         << "  resampled=" << resampled << "  rawImpacts=" << rawImpacts
								         << "  denseImpact=no"
								         << "  viRefined=no"
								         << "  offLovCheck=" << offLovCheck << "  jd=" << rep.event.jd
								         << "  missKm=" << rep.trace.distanceKm << "  sigmaLov=" << rep.trace.sigma
								         << "  p=" << std::scientific << std::setprecision(12)
								         << rep.probability.probability << std::fixed << std::setprecision(6)
								         << "  candidate=no"
								         << "  skipped=off-lov-check";
								continue;
							}
							if (checked.record) {
								rep = std::move(*checked.record);
								viRefined = checked.status == OffLovCheckStatus::Refined;
							}
						}

						if (!isReportableRisk(rep)) {
							LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] row"
							         << "  return=" << i << "  sample=" << rep.sampleIndex << "  samples=" << ret.size()
							         << "  rawSamples=" << rawSamples << "  resampled=" << resampled
							         << "  rawImpacts=" << rawImpacts
							         << "  denseImpact=" << (denseImpact ? "yes" : "no")
							         << "  viRefined=" << (viRefined ? "yes" : "no") << "  offLovCheck=" << offLovCheck
							         << "  jd=" << rep.event.jd << "  missKm=" << rep.trace.distanceKm
							         << "  distanceBodyR=" << rep.trace.distanceBodyRadii
							         << "  localMoidBodyR=" << localMoidBodyRadii(rep)
							         << "  minPossibleBodyR=" << minimumPossibleBodyRadii(rep, traceCtx.cfg.sigmaMax)
							         << "  widthBodyR=" << rep.trace.widthBodyRadii
							         << "  stretchBodyR=" << rep.trace.stretchBodyRadii
							         << "  stretchLovBodyR=" << rep.trace.stretchLovBodyRadiiPerSigma
							         << "  sigmaLov=" << rep.trace.sigma
							         << "  sigmaImpact=" << rep.probability.sigmaImpact << "  p=" << std::scientific
							         << std::setprecision(12) << rep.probability.probability << std::fixed
							         << std::setprecision(6)
							         << "  alphaDeg=" << rep.trace.alphaRad * Constants::RAD_TO_DEG
							         << "  alphaLovDeg=" << rep.trace.alphaLovRad * Constants::RAD_TO_DEG
							         << "  dD2dSigma=" << rep.trace.dDistance2DSigmaBodyRadii2 << "  candidate=no"
							         << "  skipped=no-risk";
							continue;
						}

						returnRows.push_back(makeRow(static_cast<std::size_t>(i), ret, rep, rep.covariance));

						LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] row"
						         << "  return=" << i << "  sample=" << rep.sampleIndex << "  samples=" << ret.size()
						         << "  rawSamples=" << rawSamples << "  resampled=" << resampled
						         << "  rawImpacts=" << rawImpacts << "  denseImpact=" << (denseImpact ? "yes" : "no")
						         << "  viRefined=" << (viRefined ? "yes" : "no") << "  offLovCheck=" << offLovCheck
						         << "  riskCandidate=" << (rep.riskCandidate ? "yes" : "no") << "  jd=" << rep.event.jd
						         << "  missKm=" << rep.trace.distanceKm << "  sigmaLov=" << rep.trace.sigma
						         << "  sigmaImpact=" << rep.probability.sigmaImpact << "  p=" << std::scientific
						         << std::setprecision(12) << rep.probability.probability << std::fixed
						         << std::setprecision(6) << "  candidate=" << (rep.probability.candidate ? "yes" : "no")
						         << "  widthBodyR=" << rep.trace.widthBodyRadii
						         << "  stretchBodyR=" << rep.trace.stretchBodyRadii
						         << "  stretchLovBodyR=" << rep.trace.stretchLovBodyRadiiPerSigma
						         << "  alphaDeg=" << rep.trace.alphaRad * Constants::RAD_TO_DEG
						         << "  alphaLovDeg=" << rep.trace.alphaLovRad * Constants::RAD_TO_DEG
						         << "  dD2dSigma=" << rep.trace.dDistance2DSigmaBodyRadii2;
					}

					slots[static_cast<std::size_t>(i)] = std::move(returnRows);
					refinedRows.fetch_add(slots[static_cast<std::size_t>(i)].size(), std::memory_order_relaxed);
				} catch (const std::exception& ex) {
					LOG_ERROR << "[LineOfVariation] return failed"
					          << "  return=" << i << "  firstSample=" << ret.front().sampleIndex
					          << "  lastSample=" << ret.back().sampleIndex << "  reason=" << ex.what();
					slots[static_cast<std::size_t>(i)].clear();
					activeReturns.fetch_sub(1, std::memory_order_relaxed);
					markReturnDone();
					continue;
				} catch (...) {
					LOG_ERROR << "[LineOfVariation] return failed"
					          << "  return=" << i << "  firstSample=" << ret.front().sampleIndex
					          << "  lastSample=" << ret.back().sampleIndex << "  reason=unknown-exception";
					slots[static_cast<std::size_t>(i)].clear();
					activeReturns.fetch_sub(1, std::memory_order_relaxed);
					markReturnDone();
					continue;
				}

				activeReturns.fetch_sub(1, std::memory_order_relaxed);
				markReturnDone();
			}
		}
		returnProgress.finish();

		for (auto& slot : slots)
			for (auto& row : slot)
				rows.push_back(std::move(row));

		std::sort(rows.begin(), rows.end(), [](const ImpactRiskRow& a, const ImpactRiskRow& b) {
			const double pa = a.result.probabilityImpact;
			const double pb = b.result.probabilityImpact;
			if (pa != pb)
				return pa > pb;
			return a.cluster.minDistAu < b.cluster.minDistAu;
		});

		if (cfg.maxRows >= 0 && rows.size() > static_cast<std::size_t>(cfg.maxRows)) {
			rows.resize(static_cast<std::size_t>(cfg.maxRows));
		}

		std::sort(rows.begin(), rows.end(), [](const ImpactRiskRow& a, const ImpactRiskRow& b) {
			return a.representativeEvent.jd < b.representativeEvent.jd;
		});

		return rows;
	}

}  // namespace orbital_core::hazard
