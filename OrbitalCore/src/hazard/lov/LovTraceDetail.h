#pragma once

#include "hazard/LineOfVariation.h"
#include "hazard/CloseApproachTargetPlaneDifferential.h"
#include "hazard/LovCoordinateChart.h"
#include "hazard/LovImpactProbability.h"
#include "hazard/TargetPlaneTrace.h"
#include "orbital_core/Constants.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/hazard/BPlane.h"
#include "orbital_core/hazard/CloseApproach.h"
#include "orbital_core/hazard/ImpactRisk.h"
#include "orbital_core/hazard/LineOfVariationImpact.h"

#include <Eigen/Core>

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <vector>

class IEphemerisSource;
class NBodyPropagator;
class Orbit;

namespace orbital_core::hazard::lov_detail {

	inline constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
	inline constexpr double kInvSqrtTwoPi = 0.398942280401432677939946059934381868;
	inline constexpr double kRiskProbabilityFloor = 1.0e-11;

	[[nodiscard]] inline constexpr int defaultLovWorkerCount(int maximumThreads) noexcept
	{
		return maximumThreads > 1 ? maximumThreads - 1 : 1;
	}

	using VecX = Eigen::VectorXd;
	using MatX = Eigen::MatrixXd;

	struct LovCatalogPoint {
		int sampleIndex{-1};
		double sigma{};
		VecX coordinates;
		VecX parameters;
		VecX physicalStepPerSigma;
		FitParameterCovariance coordinateCovariance;
		FitParameterCovariance physicalCovariance;
		double chiSquare{};
		bool valid{false};
	};

	struct TraceRecord {
		int sampleIndex{-1};
		VecX parameters;
		VecX stepPerSigma;
		FitParameterCovariance covariance;
		CloseApproachEvent event;
		TargetPlaneTracePoint trace;
		LovImpactProbabilityResult probability;
		bool riskCandidate{false};
		MatX dCoordsDParams;
		BPlaneBasis basis;
		double gamma{1.0};
		std::vector<double> stmCheckpointJds;
	};

	using TraceEvaluation = CloseApproachTargetPlaneTrace;

	struct LovTimingMetric {
		std::atomic<std::uint64_t> calls{0};
		std::atomic<std::uint64_t> nanoseconds{0};
	};

	struct LovTimingCounters {
		LovTimingMetric closeApproachScan;
		LovTimingMetric targetPlaneDifferential;
		LovTimingMetric minimumSearch;
		LovTimingMetric offLovRefinement;
		std::atomic<std::uint64_t> stmSegments{0};
		std::atomic<std::uint64_t> stmAcceptedSteps{0};
		std::atomic<std::uint64_t> stmRejectedSteps{0};
		std::atomic<std::uint64_t> stmRhsEvaluations{0};
		std::atomic<std::uint64_t> stmCorrectorIterations{0};
		std::atomic<std::uint64_t> stmIntegratedMicrodays{0};
	};

	struct TraceEvaluationContext {
		const LovCoordinateChart& chart;
		const FitParameterCovariance& covariance;
		const NonGravModel& nominalNonGrav;
		const Orbit& nominalOrbit;
		const NBodyPropagator& propagator;
		const LineOfVariationImpactConfig& cfg;
		IEphemerisSource& eph;
		const std::vector<LovCatalogPoint>* catalog{nullptr};
		LovTimingCounters* timings{nullptr};
	};

	enum class LovPairType {
		NoConcern,
		SimpleMinimum,
		SimpleMaximum,
		EntangledAscending,
		EntangledDescending,
		InterruptedMinimum,
		InterruptedMaximum,
		InterruptedUnknown,
		Unknown
	};

	enum class OffLovCheckStatus { NotRun, AlreadyInside, Refined, Rejected };

	struct OffLovCheckResult {
		OffLovCheckStatus status{OffLovCheckStatus::NotRun};
		std::optional<TraceRecord> record;
		const char* reason{"not-run"};
	};

	struct LovMinimumResult {
		std::optional<TraceRecord> record;
		bool converged{false};
		bool noTargetPlane{true};
		// Achilles replaces the caller's type after a derivative sign change.
		// Its INTEMIN result becomes AINTEMI, not the original INTEMIN.
		bool reclassified{false};

		[[nodiscard]] bool usable(bool requireConvergence = false) const noexcept
		{
			return record.has_value() && !noTargetPlane && (!requireConvergence || converged);
		}
	};

	using LovEncounterEvaluator =
	    std::function<std::optional<TraceRecord>(int sampleIndex, double sigma, double startJd, double endJd)>;

	// --- geometry / scaling helpers (LovTraceGeometry.cpp) ---
	[[nodiscard]] double angleDifferenceAbs(double a, double b) noexcept;
	[[nodiscard]] LovPairType classifyLovPair(const TraceRecord& left, const TraceRecord& right) noexcept;
	[[nodiscard]] double closeApproachWindowDays(double encounterRadiusAu, double vInfKmS) noexcept;
	[[nodiscard]] double minimumPossibleBodyRadii(const TraceRecord& record, double sigmaMax) noexcept;
	[[nodiscard]] double localMoidBodyRadii(const TraceRecord& record) noexcept;

	// Reference lovclosapptp/arrloadtp: use the last entered episode, reject
	// njc=0 or an out-of-window last minimum, then select its lowest minimum.
	[[nodiscard]] std::optional<std::size_t> selectReturnTargetPlaneMinimum(
	    std::span<const CloseApproachEvent> events, std::span<const std::optional<TraceEvaluation>> traces,
	    std::span<const CloseApproachPassage> passages, double windowStartJd, double windowEndJd, double epochJd,
	    double referenceJd) noexcept;
	[[nodiscard]] std::vector<double> encounterCheckpointJds(const std::vector<CloseApproachEvent>& events,
	                                                         double epochJd, double targetJd);

	// --- parameter / chart adapters (LovParameters.cpp) ---
	[[nodiscard]] VecX nominalParameters(const Orbit& orbit, const NonGravModel& ng);
	[[nodiscard]] NonGravModel nonGravFromParameters(const NonGravModel& nominal, const VecX& params);
	[[nodiscard]] LovImpactProbabilityResult probabilityForTrace(const TargetPlaneTracePoint& trace);
	[[nodiscard]] bool validConfig(const LineOfVariationImpactConfig& cfg) noexcept;
	[[nodiscard]] double nominalObservationChiSquareRaw(const TraceEvaluationContext& ctx);
	[[nodiscard]] double nominalObservationChiSquare(const TraceEvaluationContext& ctx);

	// --- reference lin_multi_gen catalog (LovCatalog.cpp) ---
	[[nodiscard]] std::vector<LovCatalogPoint> buildLovCatalog(const TraceEvaluationContext& ctx,
	                                                           const std::vector<double>& sigmas);

	[[nodiscard]] inline double catalogPredictorStep(const std::vector<double>& sigmas, std::size_t index,
	                                                 std::size_t center) noexcept
	{
		if (sigmas.empty() || index >= sigmas.size() || center >= sigmas.size() || index == center) {
			return 0.0;
		}

		// multi_gen labels point j with sigmav(j), but prop_sig advances it by
		// deltasigv(j), the interval immediately outside that label.
		if (index > center) {
			return index + 1 < sigmas.size() ? sigmas[index + 1] - sigmas[index] : sigmas[index] - sigmas[index - 1];
		}

		return index > 0 ? sigmas[index - 1] - sigmas[index] : sigmas[index] - sigmas[index + 1];
	}

	[[nodiscard]] inline double catalogSigmaQ(double chiSquare, double referenceChiSquare) noexcept
	{
		if (!std::isfinite(chiSquare) || !std::isfinite(referenceChiSquare)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const double delta = chiSquare - referenceChiSquare;
		return std::copysign(std::sqrt(std::abs(delta)), delta);
	}

	[[nodiscard]] const LovCatalogPoint* catalogAnchorAtSigma(std::span<const LovCatalogPoint> catalog,
	                                                          double sigma) noexcept;
	[[nodiscard]] std::optional<LovCatalogPoint> catalogPointAtSigma(const TraceEvaluationContext& ctx, double sigma);

	// --- per-sigma sampling / trace evaluation (LovSampling.cpp) ---
	// covariance / non-grav / propagator / eph all come from ctx; only the
	// genuinely per-call quantities are passed explicitly.
	[[nodiscard]] std::optional<TraceEvaluation> makeTraceEvaluation(
	    const TraceEvaluationContext& ctx, const CloseApproachEvent& event, const VecX& parameters,
	    const FitParameterCovariance& covariance, const VecX& stepPerSigma, double epochTdb, double sigma,
	    std::span<const double> intermediateCheckpointJds);
	[[nodiscard]] std::optional<TraceRecord> evaluateSigmaNearEncounter(int sampleIndex, double sigma,
	                                                                    double referenceStartJd, double referenceEndJd,
	                                                                    const TraceEvaluationContext& ctx);
	[[nodiscard]] std::vector<TraceRecord> evaluateSigma(const TraceEvaluationContext& ctx,
	                                                     const LovCatalogPoint& catalogPoint);
	[[nodiscard]] std::optional<TraceRecord> evaluateParametersNearEncounter(const TraceRecord& source,
	                                                                         const VecX& parameters, double referenceJd,
	                                                                         const TraceEvaluationContext& ctx,
	                                                                         double minimumWindowDays = 20.0);

	// --- bracketing / Newton / VI refinement (LovRefinement.cpp) ---
	[[nodiscard]] LovMinimumResult refineBracketedMinimum(const TraceRecord& left, const TraceRecord& right,
	                                                      const TraceEvaluationContext& ctx);
	[[nodiscard]] LovMinimumResult refineBracketedMinimum(const TraceRecord& left, const TraceRecord& right,
	                                                      const LovEncounterEvaluator& evaluate);
	[[nodiscard]] OffLovCheckResult refineVirtualImpactorRecord(const TraceRecord& representative,
	                                                            const TraceEvaluationContext& ctx);
	[[nodiscard]] LovMinimumResult refineEndpointMinimum(const TraceRecord& start, double sigmaStepLimit,
	                                                     const TraceEvaluationContext& ctx);
	[[nodiscard]] LovMinimumResult refineEndpointMinimum(const TraceRecord& start, double sigmaStepLimit,
	                                                     const LovEncounterEvaluator& evaluate);

	// --- return analysis / grouping (LovReturnBuilding.cpp) ---
	[[nodiscard]] std::vector<TraceRecord> addRefinedMinima(const std::vector<TraceRecord>& recordsIn,
	                                                        const TraceEvaluationContext& ctx);
	[[nodiscard]] std::vector<TraceRecord> addRefinedMinima(
	    const std::vector<TraceRecord>& records, double sigmaMax,
	    const std::function<double(const TraceRecord&, bool)>& endpointLimit,
	    const std::function<LovMinimumResult(const TraceRecord&, double)>& refineEndpoint,
	    const std::function<LovMinimumResult(const TraceRecord&, const TraceRecord&)>& refineBracket);
	[[nodiscard]] std::vector<std::vector<TraceRecord>> buildReturnRecords(std::vector<TraceRecord> records);
	[[nodiscard]] std::vector<TraceRecord> addDenseImpactSegments(const std::vector<TraceRecord>& records,
	                                                              const std::vector<double>& sigmas);
	[[nodiscard]] double sampleDeltaSigma(const std::vector<double>& sigmas, int sampleIndex) noexcept;
	[[nodiscard]] double endpointSigmaStepLimit(const std::vector<LovCatalogPoint>& catalog, int sampleIndex,
	                                            bool increasingHead, double fallbackStep) noexcept;
	[[nodiscard]] int countImpacts(const std::vector<TraceRecord>& records) noexcept;
	[[nodiscard]] bool isReportableRisk(const TraceRecord& record) noexcept;
	// --- row assembly (LovRowAssembly.cpp) ---
	[[nodiscard]] ImpactRiskRow makeRow(std::size_t clusterIndex, const std::vector<TraceRecord>& records,
	                                    const TraceRecord& representative, const FitParameterCovariance& covariance);
}  // namespace orbital_core::hazard::lov_detail
