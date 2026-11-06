#include "strewn_orbit/OrbitJob.h"
#include "orbital_core/hazard/SolutionPropagation.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>

namespace strewn {
	namespace {
		using Json = nlohmann::json;
		namespace hz = orbital_core::hazard;

		void keys(const Json& object, std::initializer_list<std::string_view> allowed)
		{
			if (!object.is_object())
				throw std::invalid_argument("Expected an orbit configuration object");
			for (const auto& [key, value] : object.items())
				if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
					throw std::invalid_argument("Unknown orbit configuration key: " + key);
		}

		std::uint64_t count(const Json& object, const char* name, std::uint64_t fallback)
		{
			if (!object.contains(name))
				return fallback;
			const auto& value = object.at(name);
			if (!value.is_number_integer() || (!value.is_number_unsigned() && value.get<std::int64_t>() < 0))
				throw std::invalid_argument(std::string(name) + " must be a nonnegative integer");
			return value.get<std::uint64_t>();
		}

		Json vector(const Vec3D& value) { return {value.x, value.y, value.z}; }

		std::filesystem::path allocateRun(const std::filesystem::path& root)
		{
			const auto runs = root / "runs";
			std::filesystem::create_directories(runs);
			const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(
			                      std::chrono::system_clock::now().time_since_epoch())
			                      .count();
			for (unsigned i = 0; i < 1000; ++i) {
				const auto directory = runs / ("run-" + std::to_string(time) + "-" + std::to_string(i));
				if (std::filesystem::create_directory(directory))
					return directory;
			}
			throw std::runtime_error("Cannot allocate a unique atmospheric output directory");
		}

		Json matrix(const orbital_core::FitParameterCovariance& value)
		{
			Json result = Json::array();
			for (Eigen::Index i = 0; i < value.rows(); ++i) {
				Json row = Json::array();
				for (Eigen::Index j = 0; j < value.cols(); ++j)
					row.push_back(value(i, j));
				result.push_back(std::move(row));
			}
			return result;
		}

		void write(const std::filesystem::path& file, const Json& value)
		{
			std::ofstream stream(file, std::ios::binary);
			stream.exceptions(std::ios::badbit | std::ios::failbit);
			stream << value.dump(2) << '\n';
		}

		hz::SolutionPropagationConfig propagation(const Json& input, unsigned workers)
		{
			keys(input, {"kind", "path", "clones", "seed", "altitude_km", "propagation"});
			hz::SolutionPropagationConfig config;
			config.workers = workers;
			const auto clones = count(input, "clones", 1000);
			if (clones == 0 || clones > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
				throw std::invalid_argument("Orbital clones must be in [1, INT_MAX]");
			config.sampling.numClones = static_cast<int>(clones);
			config.sampling.rngSeed = count(input, "seed", 42);
			config.sampling.altitudeKm = input.value("altitude_km", 100.0);
			if (!std::isfinite(config.sampling.altitudeKm) || config.sampling.altitudeKm <= 0)
				throw std::invalid_argument("An atmospheric entry altitude must be positive");
			const auto& p = input.at("propagation");
			keys(p, {"integrator", "horizon_days", "step_days", "planet_mask", "sb441", "encounter_radius_au",
			         "threshold_au", "sample_step_days"});
			const auto integrator = p.value("integrator", "ra15");
			if (integrator == "ra15")
				config.integrator = IntegratorKind::RA15;
			else if (integrator == "ias15")
				config.integrator = IntegratorKind::IAS15;
			else
				throw std::invalid_argument("Orbital integrator must be ra15 or ias15");
			config.forces.enabled = true;
			if (!p.contains("planet_mask") || !p.contains("sb441"))
				throw std::invalid_argument(
				    "Set planet_mask and sb441 explicitly: old solution files do not store the force model");
			const auto mask = count(p, "planet_mask", 0);
			if (mask > 0xffff)
				throw std::invalid_argument("planet_mask must fit uint16");
			config.forces.planetMask = static_cast<std::uint16_t>(mask);
			config.forces.useSb441 = p.at("sb441").get<bool>();
			config.forces.stepSize = p.value("step_days", 10.0);
			config.encounter.horizonYears = p.at("horizon_days").get<double>() / 365.25;
			config.encounter.planet = hz::Planet::Earth;
			config.encounter.thresholdAu = p.value("threshold_au", .1);
			config.encounter.encounterRadiusAu = p.value("encounter_radius_au", .1);
			config.encounter.sampleStepDay = p.value("sample_step_days", .01);
			return config;
		}

		Json orbitReport(const orbital_core::OrbitSolution& solution, const hz::SolutionPropagationConfig& config,
		                 const hz::SolutionPropagationResult& result)
		{
			const auto& mc = result.mc;
			Json events = Json::array(), crossings = Json::array();
			for (const auto& event : result.events)
				events.push_back({{"epoch_mjd_tdb", event.jd - 2400000.5}, {"distance_au", event.distAu}});
			for (const auto& point : mc.corridor)
				crossings.push_back({{"id", point.cloneIndex},
				                     {"event", point.eventIndex},
				                     {"epoch_mjd_tdb", point.epochMjdTdb},
				                     {"latitude_deg", point.latitudeDeg},
				                     {"longitude_deg", point.longitudeDeg}});
			const auto& ng = solution.nonGrav;
			return {
			    {"schema_version", 1},
			    {"object", solution.objectName},
			    {"state",
			     {{"frame", "heliocentric ECLM J2000"},
			      {"epoch_mjd_tdb", solution.nominal.epochMjd()},
			      {"position_au", vector(solution.nominal.getR())},
			      {"velocity_au_day", vector(solution.nominal.getV() * Constants::GAUSS_K)}}},
			    {"covariance", matrix(solution.covariance)},
			    {"additional_covariance_scale", 1.0},
			    {"covariance_layout", "R_AU,V_AU_day,non_grav_in_fit_mask_order"},
			    {"non_grav",
			     {{"enabled", ng.enabled()},
			      {"model", static_cast<int>(ng.model())},
			      {"fit_mask", static_cast<int>(ng.fitParameterMask())},
			      {"a1_au_day2", ng.A1()},
			      {"a2_au_day2", ng.A2()},
			      {"a3_au_day2", ng.A3()},
			      {"dt_day", ng.DT()}}},
			    {"integrator", config.integrator == IntegratorKind::RA15 ? "ra15" : "ias15"},
			    {"planet_mask", config.forces.planetMask},
			    {"sb441", config.forces.useSb441},
			    {"step_days", config.forces.stepSize},
			    {"horizon_days", config.encounter.horizonYears * 365.25},
			    {"encounter_radius_au", config.encounter.encounterRadiusAu},
			    {"threshold_au", config.encounter.thresholdAu},
			    {"sample_step_days", config.encounter.sampleStepDay},
			    {"requested", config.sampling.numClones},
			    {"seed", config.sampling.rngSeed},
			    {"workers", config.workers},
			    {"completed", mc.numConverged},
			    {"failed", mc.numFailedFit},
			    {"attempted", mc.numAttempted},
			    {"entry_probability_among_completed", mc.aggregateImpactProb},
			    {"hits_per_event", mc.hitsPerEvent},
			    {"cancelled", mc.cancelled},
			    {"entries", mc.corridor.size()},
			    {"events", events},
			    {"crossings", crossings},
			    {"altitude_km", config.sampling.altitudeKm},
			    {"boundary", "WGS84 geodetic height"},
			    {"frame_model", "OrbitalCore IAU1976/1980 UT1/GAST; no polar motion"},
			    {"conditioning",
			     "atmospheric sampling is conditional on successful Earth entries in the selected encounter windows"}};
		}
	}  // namespace

	JobResult runConfiguredJob(const std::filesystem::path& path, const RunControl& control,
	                           std::function<void(int, int)> orbitalProgress)
	{
		std::ifstream stream(path, std::ios::binary);
		if (!stream)
			throw std::runtime_error("Cannot read job configuration");
		const auto request = Json::parse(stream);
		if (request.at("input").at("kind") != "orbit_solution")
			return StrewnJob{}.run(loadJob(path), control);
		auto job = loadJob(path, EntrySamples{});
		if (!job.outputDirectory.empty() &&
		    (std::filesystem::exists(job.outputDirectory) && !std::filesystem::is_empty(job.outputDirectory)))
			throw std::invalid_argument("Output directory must be new or empty");
		const auto config = propagation(request.at("input"), job.monteCarlo.workers);
		const auto text = request.at("input").at("path").get<std::string>();
		const auto directory = std::filesystem::absolute(path).parent_path() /
		                       std::filesystem::path(std::u8string(text.begin(), text.end()));
		const auto solution = orbital_core::loadOrbitSolution(directory);
		if (control.cancel.stop_requested()) {
			JobResult result;
			result.simulation.cancelled = true;
			return result;
		}
		job.environmentDirectory = directory / "StrewnField";
		if (job.outputDirectory.empty())
			job.outputDirectory = allocateRun(job.environmentDirectory);
		const auto propagated = hz::propagateSolution(solution, config, control.cancel, [&](int done, int total) {
			if (orbitalProgress)
				orbitalProgress(done, total);
			return !control.cancel.stop_requested();
		});
		const auto report = orbitReport(solution, config, propagated);
		EntrySamples prepared;
		std::set<int> selected;
		for (const auto& point : propagated.mc.corridor) {
			if (!point.bodyFixedValid)
				throw std::runtime_error("An orbital entry has no Earth-fixed state");
			if (!selected.insert(point.cloneIndex).second)
				throw std::runtime_error("Ambiguous entry assigned to more than one encounter; narrow the horizon");
			const auto r = point.rBodyFixedKm * 1000.0, v = point.vBodyFixedKmS * 1000.0;
			prepared.values.push_back({{r.x, r.y, r.z},
			                           {v.x, v.y, v.z},
			                           point.epochMjdUtc,
			                           1.0,
			                           static_cast<std::uint64_t>(point.cloneIndex)});
		}
		const auto archiveOrbit = [&] {
			write(job.outputDirectory / "orbital_mc.json", report);
			write(job.outputDirectory / "request.json", request);
		};
		if (propagated.mc.cancelled || prepared.values.empty() || propagated.mc.numFailedFit > 0) {
			std::filesystem::create_directories(job.outputDirectory);
			archiveOrbit();
			write(job.outputDirectory / "manifest.json", {{"schema_version", 1},
			                                              {"status", propagated.mc.cancelled ? "cancelled" : "failed"},
			                                              {"stage", "orbital_entry"}});
			if (propagated.mc.cancelled) {
				JobResult result;
				result.simulation.cancelled = true;
				return result;
			}
			throw std::runtime_error(prepared.values.empty()
			                             ? "No Earth entries in the requested horizon; see orbital_mc.json"
			                             : "Orbital MC has failed clones; atmospheric sampling was not started");
		}
		job.input = std::move(prepared);
		JobResult result;
		try {
			result = StrewnJob{}.run(job, control);
		} catch (const std::exception& error) {
			if (!std::filesystem::exists(job.outputDirectory) || std::filesystem::is_empty(job.outputDirectory)) {
				std::filesystem::create_directories(job.outputDirectory);
				archiveOrbit();
				write(job.outputDirectory / "manifest.json", {{"schema_version", 1},
				                                              {"status", "failed"},
				                                              {"stage", "atmospheric_preparation"},
				                                              {"error", error.what()}});
			}
			throw;
		}
		archiveOrbit();
		if (result.simulation.cancelled && !std::filesystem::exists(job.outputDirectory / "manifest.json"))
			write(job.outputDirectory / "manifest.json",
			      {{"schema_version", 1}, {"status", "cancelled"}, {"stage", "atmospheric_preparation"}});
		return result;
	}

}  // namespace strewn
