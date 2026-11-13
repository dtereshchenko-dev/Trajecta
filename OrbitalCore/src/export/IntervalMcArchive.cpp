#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/Exporter.h"
#include "hazard/McBoundaryScan.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <array>
#include <chrono>
#include <thread>
#ifdef _WIN32
#include <Windows.h>
#endif

namespace orbital_core::hazard {
	namespace {
		using Json = nlohmann::json;

		void publishArchive(const std::filesystem::path& staging, const std::filesystem::path& destination)
		{
			std::error_code error;
			std::filesystem::rename(staging, destination, error);
#ifdef _WIN32
			constexpr std::array delays{25, 50, 100, 200, 400, 800};
			for (const int milliseconds : delays) {
				if (!error || (error.value() != ERROR_ACCESS_DENIED && error.value() != ERROR_SHARING_VIOLATION &&
				               error.value() != ERROR_LOCK_VIOLATION))
					break;
				std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
				std::filesystem::rename(staging, destination, error);
			}
#endif
			if (error)
				throw std::filesystem::filesystem_error(
				    "Cannot finalize MC archive; completed files remain in the .partial folder", staging, destination,
				    error);
		}

		Json ng(const NonGravModel& value)
		{
			return {{"model", static_cast<int>(value.model())},
			        {"enabled", value.enabled()},
			        {"fit_mask", static_cast<int>(value.fitParameterMask())},
			        {"A1", value.A1()},
			        {"A2", value.A2()},
			        {"A3", value.A3()},
			        {"DT", value.DT()}};
		}

		Json vector(const Vec3D& value) { return {value.x, value.y, value.z}; }

		Json clone(const McClone& value)
		{
			Json entry = nullptr;
			if (value.entry)
				entry = {{"epoch", value.entry->mjdTdb},
				         {"relative_position_km", vector(value.entry->relativePositionKm)},
				         {"relative_velocity_km_s", vector(value.entry->relativeVelocityKmS)},
				         {"latitude_deg", value.entry->latitudeDeg},
				         {"longitude_deg", value.entry->longitudeDeg}};
			return {{"index", value.index},
			        {"status", static_cast<int>(value.status)},
			        {"epoch", value.epochMjdTdb},
			        {"state", value.state},
			        {"non_grav", ng(value.nonGrav)},
			        {"entry_before_interval", value.entryBeforeInterval},
			        {"entry", entry},
			        {"error", value.error}};
		}
	}  // namespace

	void saveIntervalMonteCarlo(const std::filesystem::path& destination, const FitPipelineResult& source,
	                            const IntervalMcResult& result)
	{
		if (std::filesystem::exists(destination))
			throw std::runtime_error("MC archive directory already exists");
		auto directory = destination;
		directory += ".partial";
		if (!std::filesystem::create_directory(directory))
			throw std::runtime_error("MC archive directory already exists");
		const auto& c = result.config;
		const auto& forces = source.effectiveConfig.nbody;
		Json document = {{"schema", "orbital-interval-mc-4"},
		                 {"time_scale", "MJD TDB"},
		                 {"frame", "heliocentric ecliptic J2000"},
		                 {"state_units", "AU, AU/day"},
		                 {"status_codes", {"unstarted", "complete", "fit_failed", "propagation_failed", "cancelled"}},
		                 {"method", c.method == McMethod::Covariance ? "covariance" : "observations"},
		                 {"target", static_cast<int>(c.target)},
		                 {"start", c.startMjdTdb},
		                 {"end", c.endMjdTdb},
		                 {"altitude_km", c.altitudeKm},
		                 {"event_grouping", "consecutive crossing gaps; presentation only"},
		                 {"event_max_gap_days", mc::eventMaxGapDays},
		                 {"requested", c.samples},
		                 {"seed", c.seed},
		                 {"workers", c.workers},
		                 {"integrator", result.config.integrator == IntegratorKind::IAS15 ? "IAS15" : "RA15"},
		                 {"earth_frame", earthBodyFixedFrameModel()},
		                 {"body_fixed_frame", result.bodyFixedFrame},
		                 {"body_fixed_model", result.bodyFixedModel},
		                 {"forces",
		                  {{"step_days", forces.stepSize},
		                   {"planet_mask", forces.planetMask},
		                   {"sb441", forces.useSb441},
		                   {"relativity", relativityModeName(forces.forceModel.relativity)},
		                   {"sun_j2", true},
		                   {"earth_j2", forces.forceModel.earthJ2},
		                   {"earth_j3", forces.forceModel.earthJ3},
		                   {"earth_j4", forces.forceModel.earthJ4},
		                   {"earth_rotating_pole", forces.forceModel.rotatingEarthPole},
		                   {"target_naif_id", forces.targetNaifId},
		                   {"interior_mass_policy", static_cast<int>(forces.inactiveInteriorMassPolicy)}}},
		                 {"weight_profile", static_cast<int>(source.effectiveConfig.weightProfile)},
		                 {"use_weights", source.effectiveConfig.useWeights},
		                 {"accepted_indices", source.fit.keptObsIndices},
		                 {"cancelled", result.cancelled},
		                 {"completed", result.completed},
		                 {"failed_fit", result.failedFit},
		                 {"failed_propagation", result.failedPropagation},
		                 {"cancelled_clones", result.cancelledClones},
		                 {"unstarted", result.unstarted},
		                 {"entry_clones", result.entryClones},
		                 {"entries_before_interval", result.entriesBeforeInterval},
		                 {"nominal", clone(result.nominal)},
		                 {"probability_denominator", "complete trials only"},
		                 {"uncertainty_interval", "Wilson 95%, conditional on successful propagation"},
		                 {"clones", Json::array()},
		                 {"events", Json::array()},
		                 {"fit_parameter_covariance", Json::array()}};
		const auto& covariance = source.fit.uncertainty.fitParameterCovariance;
		for (Eigen::Index i = 0; i < covariance.rows(); ++i) {
			Json row = Json::array();
			for (Eigen::Index j = 0; j < covariance.cols(); ++j)
				row.push_back(covariance(i, j));
			document["fit_parameter_covariance"].push_back(std::move(row));
		}
		for (const auto& value : result.clones)
			document["clones"].push_back(clone(value));
		for (const auto& e : result.events)
			document["events"].push_back({{"start", e.startMjdTdb},
			                              {"end", e.endMjdTdb},
			                              {"entries", e.entryClones},
			                              {"probability", e.probability},
			                              {"probability_low", e.probabilityLow},
			                              {"probability_high", e.probabilityHigh}});
		std::ofstream numerical(directory / "result.json", std::ios::binary);
		numerical.exceptions(std::ios::failbit | std::ios::badbit);
		numerical << document.dump(2) << '\n';
		numerical.close();
		const auto exported =
		    Exporter::save({.objectName = "source",
		                    .fitResult = &source.fit,
		                    .observations = source.observations.all,
		                    .weightProfile = source.effectiveConfig.weightProfile,
		                    .useWeights = source.effectiveConfig.useWeights,
		                    .rootDirectory = directory},
		                   {ExportProduct::OrbitElements, ExportProduct::Observations, ExportProduct::Residuals});
		renderIntervalMonteCarloReport(directory);
		publishArchive(directory, destination);
	}
}  // namespace orbital_core::hazard
