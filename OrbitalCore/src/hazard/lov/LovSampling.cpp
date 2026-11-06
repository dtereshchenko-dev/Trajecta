#include "hazard/lov/LovTraceDetail.h"
#include <limits>
#include <utility>
#include <cstdint>

#include "NBodyPropagator.h"
#include "hazard/CloseApproachKinematics.h"
#include "hazard/EncounterPlane.h"
#include "hazard/HazardAudit.h"
#include "orbital_core/Constants.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/PlanetData.h"
#include "logger/Logger.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

namespace orbital_core::hazard::lov_detail {
	namespace {
		class ScopedTiming final {
		   public:
			explicit ScopedTiming(LovTimingMetric* metric) noexcept
			    : m_metric(metric), m_started(std::chrono::steady_clock::now())
			{
			}

			~ScopedTiming()
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

		void recordStmIntegration(const TraceEvaluationContext& ctx, const IntegrationReport& report, double epochTdb,
		                          std::span<const CloseApproachEvent> events) noexcept
		{
			if (!ctx.timings)
				return;

			ctx.timings->stmSegments.fetch_add(report.completedOutputs, std::memory_order_relaxed);
			ctx.timings->stmAcceptedSteps.fetch_add(report.acceptedSteps, std::memory_order_relaxed);
			ctx.timings->stmRejectedSteps.fetch_add(report.rejectedSteps, std::memory_order_relaxed);
			ctx.timings->stmRhsEvaluations.fetch_add(report.rhsEvaluations, std::memory_order_relaxed);
			ctx.timings->stmCorrectorIterations.fetch_add(report.correctorIterations, std::memory_order_relaxed);

			double spanDays = 0.0;
			for (const CloseApproachEvent& event : events)
				spanDays = std::max(spanDays, std::abs(event.jd - epochTdb));
			ctx.timings->stmIntegratedMicrodays.fetch_add(static_cast<std::uint64_t>(std::llround(spanDays * 1.0e6)),
			                                              std::memory_order_relaxed);
		}

		[[nodiscard]] bool isPhysicalImpact(const CloseApproachEvent& event) noexcept
		{
			const double radiusKm = planetRadiusKm(event.planet);
			return radiusKm > 0.0 && event.distAu >= 0.0 && event.distAu * kAuKm < radiusKm;
		}

		[[nodiscard]] double circularEarthVInfinityForParametersKmS(const VecX& parameters, double epochTdb) noexcept
		{
			if (parameters.size() < 6)
				return std::numeric_limits<double>::quiet_NaN();

			auto orbit = make_orbit_from_states_vectors(
			    StateVectors(Vec3D(parameters(0), parameters(1), parameters(2)),
			                 Vec3D(parameters(3), parameters(4), parameters(5)) * (1.0 / Constants::GAUSS_K)),
			    epochTdb);
			return orbit ? orbital_core::hazard::circularEarthVInfinityKmS(*orbit)
			             : std::numeric_limits<double>::quiet_NaN();
		}

		[[nodiscard]] CloseApproachConfig lovCloseApproachConfig(const CloseApproachConfig& config) noexcept
		{
			CloseApproachConfig out = config;
			out.thresholdAu = std::min(config.thresholdAu, config.encounterRadiusAu);
			return out;
		}

	}  // namespace

	std::optional<TraceEvaluation> makeTraceEvaluation(const TraceEvaluationContext& ctx,
	                                                   const CloseApproachEvent& event, const VecX& parameters,
	                                                   const FitParameterCovariance& covariance,
	                                                   const VecX& stepPerSigma, double epochTdb, double sigma,
	                                                   std::span<const double> intermediateCheckpointJds)
	{
		ScopedTiming timing(ctx.timings ? &ctx.timings->targetPlaneDifferential : nullptr);
		const CloseApproachConfig& closeApproach = ctx.cfg.closeApproach;
		CloseApproachTargetPlaneDifferentialInput input;
		input.event = &event;
		input.parameters = parameters;
		input.parameterCovariance = covariance;
		input.parameterStepPerSigma = stepPerSigma;
		input.nominalNonGrav = &ctx.nominalNonGrav;
		input.propagator = &ctx.propagator;
		input.eph = &ctx.eph;
		input.closeApproach = closeApproach;
		input.planeKind = EncounterPlaneKind::OpikTargetPlane;
		input.epochTdb = epochTdb;
		input.sigma = sigma;
		IntegrationReport integrationReport;
		input.integrationReport = &integrationReport;

		std::optional<CloseApproachTargetPlaneTrace> trace =
		    makeCloseApproachTargetPlaneTracePoint(input, intermediateCheckpointJds);
		recordStmIntegration(ctx, integrationReport, epochTdb, std::span<const CloseApproachEvent>(&event, 1));
		return trace;
	}

	namespace {
		[[nodiscard]] std::vector<std::optional<TraceEvaluation>> makeTraceEvaluations(
		    const TraceEvaluationContext& ctx, std::span<const CloseApproachEvent> events, const VecX& parameters,
		    const FitParameterCovariance& covariance, const VecX& stepPerSigma, double epochTdb, double sigma)
		{
			if (events.empty())
				return {};

			ScopedTiming timing(ctx.timings ? &ctx.timings->targetPlaneDifferential : nullptr);
			CloseApproachTargetPlaneDifferentialInput input;
			input.parameters = parameters;
			input.parameterCovariance = covariance;
			input.parameterStepPerSigma = stepPerSigma;
			input.nominalNonGrav = &ctx.nominalNonGrav;
			input.propagator = &ctx.propagator;
			input.eph = &ctx.eph;
			input.closeApproach = ctx.cfg.closeApproach;
			input.planeKind = EncounterPlaneKind::OpikTargetPlane;
			input.epochTdb = epochTdb;
			input.sigma = sigma;
			IntegrationReport integrationReport;
			input.integrationReport = &integrationReport;

			auto traces = makeCloseApproachTargetPlaneTracePoints(input, events);
			recordStmIntegration(ctx, integrationReport, epochTdb, events);
			return traces;
		}
	}  // namespace

	std::optional<TraceRecord> evaluateSigmaNearEncounter(int sampleIndex, double sigma, double referenceStartJd,
	                                                      double referenceEndJd, const TraceEvaluationContext& ctx)
	{
		audit::Scope auditScope(sampleIndex, sigma);
		const std::optional<LovCatalogPoint> point = catalogPointAtSigma(ctx, sigma);
		if (!point || !point->valid)
			return std::nullopt;

		const VecX& x = point->parameters;
		const NonGravModel ng = nonGravFromParameters(ctx.nominalNonGrav, x);

		CloseApproachConfig localCfg = lovCloseApproachConfig(ctx.cfg.closeApproach);
		const double epoch = ctx.nominalOrbit.getEpoch();
		const double vInfKmS = circularEarthVInfinityForParametersKmS(x, epoch);
		const double windowDays = closeApproachWindowDays(localCfg.encounterRadiusAu, vInfKmS);
		const double start = std::min(referenceStartJd, referenceEndJd) - windowDays;
		const double end = std::max(referenceStartJd, referenceEndJd) + windowDays;
		localCfg.lookbackYears = std::max(0.0, (epoch - start) / 365.25);
		localCfg.horizonYears = std::max(1.0e-9, (end - epoch) / 365.25);
		audit::write("encounter-input", [&](auto& log) {
			log << " startJd=" << start << " endJd=" << end << " vInfKmS=" << vInfKmS;
			for (int i = 0; i < x.size(); ++i)
				log << " x" << i << '=' << x(i);
		});

		std::vector<CloseApproachEvent> events;
		std::vector<CloseApproachPassage> passages;
		try {
			ScopedTiming timing(ctx.timings ? &ctx.timings->closeApproachScan : nullptr);
			events = findCloseApproachesFromState(Vec3D(x(0), x(1), x(2)), Vec3D(x(3), x(4), x(5)), epoch,
			                                      ctx.propagator, ng, localCfg, &ctx.eph, passages);
		} catch (const std::exception&) {
			return std::nullopt;
		} catch (...) {
			return std::nullopt;
		}

		std::vector<std::optional<TraceEvaluation>> traces =
		    makeTraceEvaluations(ctx, events, x, point->physicalCovariance, point->physicalStepPerSigma, epoch, sigma);
		const std::optional<std::size_t> selectedIndex = selectReturnTargetPlaneMinimum(
		    events, traces, passages, start, end, epoch, 0.5 * (referenceStartJd + referenceEndJd));
		audit::write("encounter-selection", [&](auto& log) {
			log << " events=" << events.size() << " passages=" << passages.size()
			    << " selected=" << (selectedIndex ? static_cast<int>(*selectedIndex) : -1);
			for (std::size_t i = 0; i < passages.size(); ++i)
				log << " passage" << i << "Group=" << passages[i].encounterIndex << " passage" << i
				    << "FirstInsideJd=" << passages[i].firstInsideJd << " passage" << i
				    << "LastInsideJd=" << passages[i].lastInsideJd;
			for (std::size_t i = 0; i < events.size(); ++i)
				log << " event" << i << "Jd=" << events[i].jd << " event" << i << "Group=" << events[i].encounterIndex
				    << " event" << i << "QAu=" << events[i].distAu << " event" << i
				    << "Tp=" << (traces[i] && traces[i]->trace.valid);
		});
		if (!selectedIndex)
			return std::nullopt;

		const CloseApproachEvent& selected = events[*selectedIndex];
		TraceEvaluation& trace = *traces[*selectedIndex];
		std::vector<double> checkpointJds = encounterCheckpointJds(events, epoch, selected.jd);

		return TraceRecord{sampleIndex,
		                   x,
		                   point->physicalStepPerSigma,
		                   point->physicalCovariance,
		                   trace.event,
		                   trace.trace,
		                   probabilityForTrace(trace.trace),
		                   true,
		                   std::move(trace.dCoordsDParams),
		                   trace.basis,
		                   trace.gamma,
		                   std::move(checkpointJds)};
	}

	std::vector<TraceRecord> evaluateSigma(const TraceEvaluationContext& ctx, const LovCatalogPoint& catalogPoint)
	{
		audit::Scope auditScope(catalogPoint.sampleIndex, catalogPoint.sigma);
		std::vector<TraceRecord> out;
		if (!catalogPoint.valid)
			return out;

		const VecX& x = catalogPoint.parameters;
		audit::write("catalog", [&](auto& log) {
			log << " epochJd=" << ctx.nominalOrbit.getEpoch() << " chiSquare=" << catalogPoint.chiSquare;
			for (int i = 0; i < x.size(); ++i) {
				log << " x" << i << '=' << x(i) << " step" << i << '=' << catalogPoint.physicalStepPerSigma(i);
				for (int j = 0; j < x.size(); ++j)
					log << " c" << i << j << '=' << catalogPoint.physicalCovariance(i, j);
			}
			for (int i = 0; i < catalogPoint.coordinates.size(); ++i)
				log << " equ" << i << '=' << catalogPoint.coordinates(i);
		});
		const NonGravModel ng = nonGravFromParameters(ctx.nominalNonGrav, x);

		const CloseApproachConfig localCfg = lovCloseApproachConfig(ctx.cfg.closeApproach);
		std::vector<CloseApproachEvent> events;
		try {
			ScopedTiming timing(ctx.timings ? &ctx.timings->closeApproachScan : nullptr);
			events = findCloseApproachesFromState(Vec3D(x(0), x(1), x(2)), Vec3D(x(3), x(4), x(5)),
			                                      ctx.nominalOrbit.getEpoch(), ctx.propagator, ng, localCfg, &ctx.eph);
		} catch (const std::exception& error) {
			audit::write("scan-exception", [&](auto& log) { log << " error=" << std::quoted(error.what()); });
			return out;
		} catch (...) {
			audit::write("scan-exception", [](auto& log) { log << " error=unknown"; });
			return out;
		}

		audit::write("scan-complete", [&](auto& log) { log << " events=" << events.size(); });
		if (events.empty() && audit::current) {
			try {
				const auto surface = ctx.propagator.propagateUntilSurface(
				    Vec3D(x(0), x(1), x(2)), Vec3D(x(3), x(4), x(5)), ctx.nominalOrbit.getEpoch(),
				    ctx.nominalOrbit.getEpoch() + localCfg.horizonYears * 365.25, toEphemerisPlanet(localCfg.planet),
				    planetRadiusKm(localCfg.planet), ctx.eph, ng, NBodyPropagator::SurfaceShape::Sphere);
				audit::write("surface-replay",
				             [&](auto& log) { log << " hit=" << surface.hitSurface << " jd=" << surface.jd; });
			} catch (const std::exception& error) {
				audit::write("surface-replay-failed",
				             [&](auto& log) { log << " error=" << std::quoted(error.what()); });
			}
		}

		std::size_t eventCount = events.size();
		for (std::size_t i = 0; i < events.size(); ++i) {
			if (isPhysicalImpact(events[i])) {
				eventCount = i + 1;
				break;
			}
		}

		const std::span<const CloseApproachEvent> selectedEvents(events.data(), eventCount);
		std::vector<std::optional<TraceEvaluation>> traces =
		    makeTraceEvaluations(ctx, selectedEvents, x, catalogPoint.physicalCovariance,
		                         catalogPoint.physicalStepPerSigma, ctx.nominalOrbit.getEpoch(), catalogPoint.sigma);

		out.reserve(eventCount);
		for (std::size_t i = 0; i < eventCount; ++i) {
			const CloseApproachEvent& event = events[i];
			const bool impact = isPhysicalImpact(event);
			if (!traces[i]) {
				audit::write("tp-failed", [&](auto& log) { log << " jd=" << event.jd; });
				if (impact)
					break;
				continue;
			}

			TraceEvaluation& trace = *traces[i];
			audit::write("tp-record", [&](auto& log) {
				const auto& tp = trace.trace;
				log << " jd=" << tp.jd << " xiKm=" << tp.coordsKm.x() << " zetaKm=" << tp.coordsKm.y()
				    << " stretchKm=" << tp.stretchKm << " widthKm=" << tp.widthKm << " gamma=" << trace.gamma
				    << " vInfKmS=" << trace.vInfKmS << " impact=" << tp.impact << " valid=" << tp.valid;
			});

			out.push_back({catalogPoint.sampleIndex, x, catalogPoint.physicalStepPerSigma,
			               catalogPoint.physicalCovariance, trace.event, trace.trace, probabilityForTrace(trace.trace),
			               false, std::move(trace.dCoordsDParams), trace.basis, trace.gamma,
			               encounterCheckpointJds(events, ctx.nominalOrbit.getEpoch(), event.jd)});

			if (impact)
				break;
		}
		audit::write("va-complete", [&](auto& log) { log << " records=" << out.size(); });
		return out;
	}

	std::optional<TraceRecord> evaluateParametersNearEncounter(const TraceRecord& source, const VecX& parameters,
	                                                           double referenceJd, const TraceEvaluationContext& ctx,
	                                                           double minimumWindowDays)
	{
		if (parameters.size() < 6 || source.stepPerSigma.size() != parameters.size()) {
			return std::nullopt;
		}

		const NonGravModel ng = nonGravFromParameters(ctx.nominalNonGrav, parameters);

		CloseApproachConfig localCfg = lovCloseApproachConfig(ctx.cfg.closeApproach);
		const double epoch = ctx.nominalOrbit.getEpoch();
		const double vInfKmS = circularEarthVInfinityForParametersKmS(parameters, epoch);
		const double windowDays =
		    std::max(minimumWindowDays, closeApproachWindowDays(localCfg.encounterRadiusAu, vInfKmS));
		const double start = referenceJd - windowDays;
		const double end = referenceJd + windowDays;
		localCfg.lookbackYears = std::max(0.0, (epoch - start) / 365.25);
		localCfg.horizonYears = std::max(1.0e-9, (end - epoch) / 365.25);

		std::vector<CloseApproachEvent> events;
		try {
			ScopedTiming timing(ctx.timings ? &ctx.timings->closeApproachScan : nullptr);
			events = findCloseApproachesFromState(Vec3D(parameters(0), parameters(1), parameters(2)),
			                                      Vec3D(parameters(3), parameters(4), parameters(5)), epoch,
			                                      ctx.propagator, ng, localCfg, &ctx.eph);
		} catch (...) {
			return std::nullopt;
		}

		auto selected = events.end();
		for (auto it = events.begin(); it != events.end(); ++it) {
			if (it->jd < start || it->jd > end)
				continue;
			if (selected == events.end() || it->distAu < selected->distAu) {
				selected = it;
			}
		}
		if (selected == events.end())
			return std::nullopt;

		std::vector<double> checkpointJds = encounterCheckpointJds(events, epoch, selected->jd);
		std::optional<TraceEvaluation> trace =
		    makeTraceEvaluation(ctx, *selected, parameters, source.covariance, source.stepPerSigma, epoch,
		                        source.trace.sigma, checkpointJds);
		if (!trace)
			return std::nullopt;

		TraceRecord out = source;
		out.parameters = parameters;
		out.event = trace->event;
		out.trace = trace->trace;
		out.covariance = source.covariance;
		out.probability = probabilityForTrace(trace->trace);
		out.riskCandidate = true;
		out.dCoordsDParams = trace->dCoordsDParams;
		out.basis = trace->basis;
		out.gamma = trace->gamma;
		out.stmCheckpointJds = std::move(checkpointJds);
		return out;
	}

}  // namespace orbital_core::hazard::lov_detail
