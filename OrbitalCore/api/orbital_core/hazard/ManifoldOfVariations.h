#pragma once

#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include <limits>
#include <span>

namespace orbital_core::hazard {

	enum class MovSampling { Automatic, Grid, Cobweb };
	enum class MovNodeStatus {
		Unstarted,
		OutsideRegion,
		FitFailed,
		Rejected,
		EarthBound,
		Complete,
		PropagationFailed,
		Cancelled
	};

	struct MovConfig {
		MovSampling sampling{MovSampling::Automatic};
		Planet target{Planet::Earth};
		double startMjdTdb{}, endMjdTdb{}, altitudeKm{100};
		int preliminaryNodes{40}, rangeNodes{80}, rateNodes{80};
		double sigmaMax{5}, absoluteMagnitudeMax{34.5};
		unsigned workers{15};
	};

	struct MovNode {
		int row{}, column{};
		MovNodeStatus status{MovNodeStatus::Unstarted};
		std::array<double, 6> attributable{}, state{};
		double chiSquare{std::numeric_limits<double>::quiet_NaN()};
		double chi{std::numeric_limits<double>::quiet_NaN()};
		double jacobianMov{}, jacobianSampling{}, weight{};
		std::optional<double> absoluteMagnitude;
		bool outerBoundary{}, entryBeforeInterval{};
		std::optional<McBoundaryEntry> entry;
		std::string error;
	};

	struct MovResult {
		MovConfig config;
		MovSampling sampling{};
		double epochMjdTdb{}, curvatureSnr{};
		std::array<double, 4> attributable{};
		std::string observerCode, bodyFixedFrame, bodyFixedModel;
		std::vector<MovNode> preliminary, nodes;
		bool preliminaryLogRange{}, logRange{};
		// Weighted quadrature on the sampled chi-limited MOV, not MC frequencies.
		// Propagation failures retain their weight: [lower,upper] brackets that
		// known missing mass, but does not bound unresolved fits or grid truncation.
		std::optional<double> probabilityLower, probabilityUpper;
		double outerBoundaryWeight{};
		bool cancelled{};
		std::vector<std::string> warnings;
	};

	// Ground-optical, gravity-only RA15 systematic ranging, for arcs and prediction
	// horizons no longer than 30 days. Photometry is required by the H selection.
	// A nominal is optional; only the
	// automatic/cobweb branch uses its covariance, with geodesic-curvature SNR > 3.
	// Input observations are already prepared/debiased, JD-sorted, and selected by
	// the caller; their mask stays fixed at every node. No IOD or clipping is run.
	// Epochs are TDB. ATT angles/rates are radians and radians/day, range/rate AU
	// and AU/day. The observer is evaluated at emission epoch + range/light speed,
	// as defined by the ATT chart. All nonlinear fits reuse the shared LM model.
	// MOV tangents use the formal weighted normal inverse, not the scaled posterior
	// covariance. Earth-bound nodes keep probability weight but are not propagated.
	ORBITAL_CORE_API MovResult computeManifoldOfVariations(std::span<const ObsNode> observations,
	                                                       const FitConfig& forcesAndWeights, const MovConfig& config,
	                                                       const FitResult* nominal = nullptr,
	                                                       std::stop_token cancel = {},
	                                                       std::function<void(int, int)> progress = {});

	// Writes a new numerical archive, including all nodes and effective weights;
	// rendering is separate so an HTML failure cannot discard the calculation.
	ORBITAL_CORE_API void saveManifoldOfVariations(const MovResult& result, std::span<const ObsNode> observations,
	                                               const FitConfig& settings, const std::filesystem::path& directory,
	                                               const std::string& designation);
	ORBITAL_CORE_API void renderManifoldOfVariationsReport(const std::filesystem::path& directory);

}  // namespace orbital_core::hazard
