#include "strewn_orbit/OrbitJob.h"
#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/TimeScale.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
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

		struct OrbitalPropagation {
			orbital_core::NBodySettings forces;
			hz::IntervalMcConfig sampling;
			double horizonDays{};
		};

		OrbitalPropagation propagation(const Json& input, unsigned workers)
		{
			keys(input, {"kind", "path", "clones", "seed", "altitude_km", "propagation"});
			OrbitalPropagation config;
			config.sampling.workers = workers;
			const auto clones = count(input, "clones", 1000);
			if (clones == 0 || clones > 1000000)
				throw std::invalid_argument("Orbital clones must be in [1, 1000000]");
			config.sampling.samples = static_cast<int>(clones);
			config.sampling.seed = count(input, "seed", 42);
			config.sampling.altitudeKm = input.value("altitude_km", 100.0);
			if (!std::isfinite(config.sampling.altitudeKm) || config.sampling.altitudeKm <= 0)
				throw std::invalid_argument("An atmospheric entry altitude must be positive");
			const auto& p = input.at("propagation");
			keys(p, {"integrator", "horizon_days", "step_days", "planet_mask", "sb441"});
			const auto integrator = p.value("integrator", "ra15");
			if (integrator == "ra15")
				config.sampling.integrator = IntegratorKind::RA15;
			else if (integrator == "ias15")
				config.sampling.integrator = IntegratorKind::IAS15;
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
			config.horizonDays = p.at("horizon_days").get<double>();
			if (!std::isfinite(config.horizonDays) || config.horizonDays <= 0)
				throw std::invalid_argument("horizon_days must be finite and positive");
			return config;
		}

		Json orbitReport(const orbital_core::OrbitSolution& solution, const OrbitalPropagation& config,
		                 const hz::IntervalMcResult& mc)
		{
			Json events = Json::array(), crossings = Json::array();
			for (const auto& event : mc.events)
				events.push_back({{"start_mjd_tdb", event.startMjdTdb},
				                  {"end_mjd_tdb", event.endMjdTdb},
				                  {"entries", event.entryClones},
				                  {"probability", event.probability}});
			for (const auto& clone : mc.clones)
				if (clone.status == hz::McCloneStatus::Complete && clone.entry)
					crossings.push_back({{"id", clone.index},
					                     {"epoch_mjd_tdb", clone.entry->mjdTdb},
					                     {"latitude_deg", clone.entry->latitudeDeg},
					                     {"longitude_deg", clone.entry->longitudeDeg}});
			const auto& ng = solution.nonGrav;
			return {{"schema_version", 2},
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
			        {"integrator", config.sampling.integrator == IntegratorKind::RA15 ? "ra15" : "ias15"},
			        {"planet_mask", config.forces.planetMask},
			        {"sb441", config.forces.useSb441},
			        {"step_days", config.forces.stepSize},
			        {"horizon_days", config.horizonDays},
			        {"requested", config.sampling.samples},
			        {"seed", config.sampling.seed},
			        {"workers", config.sampling.workers},
			        {"completed", mc.completed},
			        {"failed", mc.failedFit + mc.failedPropagation},
			        {"entry_probability_among_completed",
			         mc.completed ? Json(double(mc.entryClones) / mc.completed) : Json(nullptr)},
			        {"cancelled", mc.cancelled},
			        {"entries", mc.entryClones},
			        {"events", events},
			        {"crossings", crossings},
			        {"altitude_km", config.sampling.altitudeKm},
			        {"boundary", "WGS84 geodetic height"},
			        {"frame_model", "OrbitalCore IAU1976/1980 UT1/GAST; no polar motion"},
			        {"conditioning",
			         "atmospheric sampling is conditional on successful Earth entries in the requested interval"}};
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
		auto config = propagation(request.at("input"), job.monteCarlo.workers);
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
		config.sampling.startMjdTdb = solution.nominal.epochMjd();
		config.sampling.endMjdTdb = config.sampling.startMjdTdb + config.horizonDays;
		const auto propagated = hz::runIntervalMonteCarlo(solution, config.forces, config.sampling, control.cancel,
		                                                  [&](int done, int total) {
			                                                  if (orbitalProgress)
				                                                  orbitalProgress(done, total);
		                                                  });
		const auto report = orbitReport(solution, config, propagated);
		EntrySamples prepared;
		std::set<int> selected;
		for (const auto& clone : propagated.clones) {
			if (clone.status != hz::McCloneStatus::Complete || !clone.entry)
				continue;
			const auto& entry = *clone.entry;
			if (!selected.insert(clone.index).second)
				throw std::runtime_error("Duplicate orbital entry");
			const auto r = hz::icrfToBodyFixedMjd(hz::Planet::Earth, entry.relativePositionKm, entry.mjdTdb) * 1000.0;
			const auto v = hz::bodyFixedVelocityMjd(hz::Planet::Earth, entry.relativePositionKm,
			                                        entry.relativeVelocityKmS, entry.mjdTdb) *
			               1000.0;
			if (!std::isfinite(r.dot(r)) || !std::isfinite(v.dot(v)))
				throw std::runtime_error("An orbital entry has no valid Earth-fixed state");
			prepared.values.push_back({{r.x, r.y, r.z},
			                           {v.x, v.y, v.z},
			                           orbital_core::time::fromMjdTdb(entry.mjdTdb, orbital_core::time::Scale::UTC),
			                           1.0,
			                           static_cast<std::uint64_t>(clone.index)});
		}
		const auto archiveOrbit = [&] {
			write(job.outputDirectory / "orbital_mc.json", report);
			write(job.outputDirectory / "request.json", request);
		};
		if (propagated.cancelled || prepared.values.empty() || propagated.failedFit > 0 ||
		    propagated.failedPropagation > 0) {
			std::filesystem::create_directories(job.outputDirectory);
			archiveOrbit();
			write(job.outputDirectory / "manifest.json", {{"schema_version", 1},
			                                              {"status", propagated.cancelled ? "cancelled" : "failed"},
			                                              {"stage", "orbital_entry"}});
			if (propagated.cancelled) {
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
