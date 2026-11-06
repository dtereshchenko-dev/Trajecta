#include "orbital_core/hazard/MonteCarloArchive.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "orbital_core/hazard/PlanetData.h"
#include "NBodyPropagator.h"
#include <Eigen/Eigenvalues>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace orbital_core::hazard {
	namespace {
		using Json = nlohmann::json;

		Json vector(const Vec3D& v) { return {v.x, v.y, v.z}; }

		bool finite(const Vec3D& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

		template <class Matrix>
		Json matrix(const Matrix& value)
		{
			Json rows = Json::array();
			for (Eigen::Index i = 0; i < value.rows(); ++i) {
				Json row = Json::array();
				for (Eigen::Index j = 0; j < value.cols(); ++j)
					row.push_back(value(i, j));
				rows.push_back(std::move(row));
			}
			return rows;
		}

		void write(const std::filesystem::path& file, const Json& value)
		{
			std::ofstream stream(file, std::ios::binary);
			stream.exceptions(std::ios::failbit | std::ios::badbit);
			stream << value.dump(2) << '\n';
		}

		std::filesystem::path allocate(const std::filesystem::path& solutionDirectory)
		{
			if (solutionDirectory.empty() || !std::filesystem::is_directory(solutionDirectory))
				throw std::invalid_argument("MC archive requires the saved solution directory");
			const auto root = std::filesystem::absolute(solutionDirectory) / "MC";
			std::filesystem::create_directories(root);
			const auto now = std::chrono::duration_cast<std::chrono::milliseconds>(
			                     std::chrono::system_clock::now().time_since_epoch())
			                     .count();
			for (unsigned i = 0; i < 1000; ++i) {
				const auto path = root / ("run-" + std::to_string(now) + "-" + std::to_string(i));
				if (std::filesystem::create_directory(path))
					return path;
			}
			throw std::runtime_error("Cannot allocate MC archive directory");
		}

		bool validCovariance(const FitParameterCovariance& covariance)
		{
			if (!covariance.allFinite() || covariance.norm() == 0 || covariance.rows() != covariance.cols())
				return false;
			if ((covariance - covariance.transpose()).norm() > 1e-10 * covariance.norm())
				return false;
			Eigen::MatrixXd correlation = covariance;
			for (Eigen::Index i = 0; i < covariance.rows(); ++i) {
				if (covariance(i, i) <= 0)
					return false;
				for (Eigen::Index j = 0; j < covariance.cols(); ++j) {
					if (covariance(j, j) <= 0)
						return false;
					correlation(i, j) /= std::sqrt(covariance(i, i)) * std::sqrt(covariance(j, j));
				}
			}
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(correlation);
			return eigen.info() == Eigen::Success && eigen.eigenvalues().minCoeff() >= -1e-10;
		}

	}  // namespace

	ExportReport saveMonteCarloArchive(const std::filesystem::path& solutionDirectory,
	                                   const MonteCarloArchiveInput& input, const NBodyPropagator& propagator,
	                                   std::span<const CloseApproachEvent> events, const ObsMcResult& result)
	{
		const auto& s = input.solution;
		const auto& ng = s.nonGrav;
		const auto& c = input.sampling;
		const bool parameters = input.kind == CovarianceSampling::Parameters;
		const int dimension = 6 + (parameters ? ng.numFitParams() : 0);
		if (c.numClones <= 0 || !std::isfinite(c.altitudeKm) || c.altitudeKm < 0 ||
		    !std::isfinite(input.covarianceObjectiveScale) || input.covarianceObjectiveScale <= 0 ||
		    !finite(s.nominal.getR()) || !finite(s.nominal.getV()) || !std::isfinite(s.nominal.epochMjd()))
			throw std::invalid_argument("Invalid MC archive input");
		Json order = {"x_au", "y_au", "z_au", "vx_au_day", "vy_au_day", "vz_au_day"};
		if (parameters) {
			const std::array<const char*, 4> names{"a1_au_day2", "a2_au_day2", "a3_au_day2", "dt_day"};
			for (int i = 0; i < ng.numFitParams(); ++i)
				order.push_back(names[static_cast<int>(ng.fitParameterAt(i))]);
		}
		Json reasons = Json::array();
		const auto reject = [&](const char* text) {
			if (std::find(reasons.begin(), reasons.end(), text) == reasons.end())
				reasons.push_back(text);
		};
		const bool covarianceValid =
		    s.covariance.rows() == dimension && s.covariance.cols() == dimension && validCovariance(s.covariance);
		if (!covarianceValid)
			reject("invalid_input_covariance");
		if (!parameters && ng.numFitParams() > 0)
			reject("fitted_non_grav_covariance_not_sampled");
		if (result.cancelled)
			reject("cancelled");
		if (result.numAttempted != c.numClones || result.numConverged != c.numClones)
			reject("incomplete_clones");
		if (result.numFailedFit != 0)
			reject("failed_clones");
		if (result.numAnyImpactClones != c.numClones)
			reject("not_all_clones_crossed");
		if (c.altitudeKm != 100 || result.altitudeKm != c.altitudeKm)
			reject("requires_100_km_entry_boundary");
		if (result.hitsPerEvent.size() != events.size())
			reject("inconsistent_event_counts");

		Json selectedEvents = Json::array();
		for (const auto& e : events)
			selectedEvents.push_back({{"planet", static_cast<int>(e.planet)},
			                          {"epoch_mjd_tdb", e.jd - 2400000.5},
			                          {"distance_au", e.distAu},
			                          {"position_au", vector(e.rTarget)},
			                          {"velocity_au_day", vector(e.vTarget)},
			                          {"planet_position_au", vector(e.rPlanet)},
			                          {"planet_velocity_au_day", vector(e.vPlanet)}});
		Json request = {{"schema_version", 1},
		                {"sampling_kind", parameters ? "parameter_covariance" : "state_covariance"},
		                {"object", s.objectName},
		                {"nominal",
		                 {{"frame", "heliocentric ECLM J2000"},
		                  {"epoch_mjd_tdb", s.nominal.epochMjd()},
		                  {"position_au", vector(s.nominal.getR())},
		                  {"velocity_au_day", vector(s.nominal.getV() * Constants::GAUSS_K)}}},
		                {"covariance",
		                 {{"matrix", matrix(s.covariance)},
		                  {"parameter_order", order},
		                  {"objective_scale_already_applied", input.covarianceObjectiveScale},
		                  {"additional_scale", 1.0},
		                  {"valid", covarianceValid}}},
		                {"non_grav",
		                 {{"enabled", ng.enabled()},
		                  {"model", static_cast<int>(ng.model())},
		                  {"fit_mask", static_cast<int>(ng.fitParameterMask())},
		                  {"a1_au_day2", ng.A1()},
		                  {"a2_au_day2", ng.A2()},
		                  {"a3_au_day2", ng.A3()},
		                  {"dt_day", ng.DT()}}},
		                {"monte_carlo",
		                 {{"clones", c.numClones},
		                  {"seed", c.rngSeed},
		                  {"altitude_km", c.altitudeKm},
		                  {"rng", "SplitMix64 keyed Box-Muller; clone ID and parameter index"},
		                  {"noise_sigma_mult", c.noiseSigmaMult},
		                  {"noise_sigma_mult_applies", false}}},
		                {"propagation",
		                 {{"integrator", propagator.integratorName()},
		                  {"accuracy", propagator.fixedAccuracyParameter()},
		                  {"step_days", propagator.stepSize()},
		                  {"planet_mask", propagator.planetMask()},
		                  {"sb441", propagator.useSb441()},
		                  {"target_naif_id", propagator.targetNaifId()},
		                  {"encounter_points", propagator.encounterStepPointCount()},
		                  {"earth_encounter_radius_au", propagator.earthEncounterRadiusAu()},
		                  {"inactive_interior_mass_policy", static_cast<int>(propagator.inactiveInteriorMassPolicy())},
		                  {"event_window_days", 1.0},
		                  {"events", selectedEvents}}},
		                {"external_ephemerides_embedded", false}};
		request["physical"] = {{"absolute_magnitude", s.absoluteMagnitude ? Json(*s.absoluteMagnitude) : Json(nullptr)},
		                       {"albedo", input.albedo ? Json(*input.albedo) : Json(nullptr)}};
#ifdef _OPENMP
		request["monte_carlo"]["openmp_max_threads"] = omp_get_max_threads();
#else
		request["monte_carlo"]["openmp_max_threads"] = 1;
#endif

		Json entries = Json::array(), crossings = Json::array();
		std::set<int> cloneIds, eventIds;
		std::vector<int> counted(events.size());
		for (const auto& p : result.corridor) {
			crossings.push_back({{"clone", p.cloneIndex},
			                     {"event", p.eventIndex},
			                     {"planet", static_cast<int>(p.planet)},
			                     {"epoch_mjd_tdb", p.epochMjdTdb},
			                     {"epoch_mjd_utc", p.epochMjdUtc},
			                     {"body_fixed_valid", p.bodyFixedValid},
			                     {"position_body_fixed_km", vector(p.rBodyFixedKm)},
			                     {"velocity_body_fixed_km_s", vector(p.vBodyFixedKmS)},
			                     {"latitude_deg", p.latitudeDeg},
			                     {"longitude_deg", p.longitudeDeg},
			                     {"diagnostic_miss_km", p.missKm},
			                     {"relative_speed_km_s", p.velocityKmS}});
			if (p.cloneIndex < 0 || p.cloneIndex >= c.numClones || p.eventIndex < 0 ||
			    static_cast<std::size_t>(p.eventIndex) >= events.size()) {
				reject("invalid_clone_or_event_id");
				continue;
			}
			++counted[p.eventIndex];
			eventIds.insert(p.eventIndex);
			if (!cloneIds.insert(p.cloneIndex).second)
				reject("multiple_crossings_per_clone");
			if (p.planet != Planet::Earth || events[p.eventIndex].planet != Planet::Earth) {
				reject("non_earth_crossing");
				continue;
			}
			if (!p.bodyFixedValid || !finite(p.rBodyFixedKm) || !finite(p.vBodyFixedKmS) ||
			    !std::isfinite(p.epochMjdUtc) || !std::isfinite(p.epochMjdTdb)) {
				reject("invalid_entry_state");
				continue;
			}
			const auto geo = geodeticOnReferenceEllipsoidKm(p.rBodyFixedKm, Constants::EARTH_EQ_RADIUS_KM,
			                                                Constants::EARTH_FLATTENING);
			const double height = geo.altitudeKm;
			if (c.altitudeKm > 0 && (!std::isfinite(height) || std::abs(height - c.altitudeKm) > .001))
				reject("entry_height_mismatch");
			const Vec3D up{std::cos(geo.latitudeRad) * std::cos(geo.longitudeRad),
			               std::cos(geo.latitudeRad) * std::sin(geo.longitudeRad), std::sin(geo.latitudeRad)};
			if (!geo.valid || up.dot(p.vBodyFixedKmS) >= 0)
				reject("entry_not_inbound");
			const auto r = p.rBodyFixedKm * 1000, v = p.vBodyFixedKmS * 1000;
			const auto id = static_cast<std::uint64_t>(p.cloneIndex) * events.size() + p.eventIndex;
			entries.push_back({{"id", id},
			                   {"epoch_mjd_utc", p.epochMjdUtc},
			                   {"position_ecef_m", vector(r)},
			                   {"velocity_ecef_mps", vector(v)}});
		}
		if (counted != result.hitsPerEvent)
			reject("inconsistent_event_counts");
		if (cloneIds.size() != static_cast<std::size_t>(result.numAnyImpactClones))
			reject("inconsistent_clone_count");
		if (eventIds.size() != 1)
			reject("requires_one_encounter");
		if (entries.size() != static_cast<std::size_t>(c.numClones))
			reject("incomplete_entry_states");
		const bool entryReady = reasons.empty();
		const bool sizeKnown = s.absoluteMagnitude && std::isfinite(*s.absoluteMagnitude) && input.albedo &&
		                       std::isfinite(*input.albedo) && *input.albedo > 0 && *input.albedo <= 1;
		if (!sizeKnown)
			reject("missing_H_or_albedo");
		const bool ready = reasons.empty();

		ExportReport report;
		report.directory = allocate(solutionDirectory);
		const auto save = [&](const char* name, const Json& value) {
			write(report.directory / name, value);
			report.files.push_back(report.directory / name);
		};

		save("request.json", request);
		save("result.json", {{"schema_version", 1},
		                     {"attempted", result.numAttempted},
		                     {"completed", result.numConverged},
		                     {"failed", result.numFailedFit},
		                     {"cancelled", result.cancelled},
		                     {"any_crossing_clones", result.numAnyImpactClones},
		                     {"hits_per_event", result.hitsPerEvent},
		                     {"empirical_crossing_probability", result.aggregateImpactProb},
		                     {"altitude_km", result.altitudeKm},
		                     {"crossings", crossings}});
		save("entries.json",
		     {{"schema_version", 1},
		      {"frame", "ECEF"},
		      {"altitude_km", result.altitudeKm},
		      {"boundary", result.altitudeKm > 0 ? "WGS84 geodetic height" : "legacy sphere"},
		      {"frame_model", result.earthFrameModel.empty() ? "unknown (not recorded)" : result.earthFrameModel},
		      {"entries", entries}});
		save("handoff.json", {{"schema_version", 1},
		                      {"kind", "atmospheric_input"},
		                      {"ready", ready},
		                      {"automatic_launch", false},
		                      {"input", {{"kind", "state_samples"}, {"path", "entries.json"}}},
		                      {"physical", request["physical"]},
		                      {"source_request", "request.json"},
		                      {"blocked_reasons", reasons}});
		save(
		    "manifest.json",
		    {{"schema_version", 1},
		     {"status", result.cancelled
		                    ? "cancelled"
		                    : (result.numFailedFit > 0 || result.numConverged != c.numClones ? "partial" : "complete")},
		     {"entry_states_ready", entryReady},
		     {"atmospheric_input_ready", ready},
		     {"blocked_reasons", reasons},
		     {"automatic_launch", false}});
		return report;
	}

}  // namespace orbital_core::hazard
