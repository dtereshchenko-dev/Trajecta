#include "JobFiles.h"
#include "strewn_core/AtmosphericInput.h"
#include "AutoEnvironment.h"
#include "Convergence.h"
#include "../output/FieldOutputInternal.h"
#include <limits>
#include <optional>
#include <stdexcept>

namespace strewn {
	using namespace jobfiles;

	namespace {

		std::filesystem::path resolve(const std::filesystem::path& base, const Json& value)
		{
			const auto text = value.get<std::string>();
			return std::filesystem::absolute(base / std::filesystem::path(std::u8string(text.begin(), text.end())))
			    .lexically_normal();
		}

		Range range(const Json& object, const char* key, Range fallback)
		{
			if (!object.contains(key))
				return fallback;
			const auto& value = object.at(key);
			if (value.is_number()) {
				const auto scalar = value.get<double>();
				return {scalar, scalar};
			}
			const auto values = numbers<2>(value);
			return {values[0], values[1]};
		}

		template <class T>
		std::array<T, 2> counts(const Json& value)
		{
			if (!value.is_array() || value.size() != 2)
				throw std::invalid_argument("Grouping counts require exactly two integers");
			std::array<T, 2> result;
			for (std::size_t i = 0; i < result.size(); ++i) {
				if (!value[i].is_number_integer() ||
				    (!value[i].is_number_unsigned() && value[i].get<std::int64_t>() < 0) ||
				    value[i].get<std::uint64_t>() > std::numeric_limits<T>::max())
					throw std::invalid_argument("Grouping count is outside the unsigned integer range");
				result[i] = value[i].get<T>();
			}
			return result;
		}

	}  // namespace

	static JobSpec loadJobImpl(const std::filesystem::path& path, std::optional<EntrySamples> prepared)
	{
		const auto value = readJson(path);
		keys(value, {"schema_version", "input", "environment", "physics", "integration", "monte_carlo",
		             "output_directory", "environment_directory", "maps"});
		if (value.at("schema_version") != 1 && value.at("schema_version") != 2)
			throw std::invalid_argument("Unsupported job schema");
		JobSpec result;
		const bool legacy = value.at("schema_version") == 1;
		if (legacy)
			result.physics.retainMassRemainder = false;
		const auto base = std::filesystem::absolute(path).parent_path();
		if (value.contains("environment_directory"))
			result.environmentDirectory = resolve(base, value.at("environment_directory"));
		if (value.contains("maps") && !value.at("maps").is_null())
			result.maps = output::readOptions(value.at("maps"), base);
		const auto& input = value.at("input");
		const auto kind = input.at("kind").get<std::string>();
		if (value.contains("output_directory"))
			result.outputDirectory = resolve(base, value.at("output_directory"));
		else if (!prepared || kind != "orbit_solution")
			throw std::invalid_argument("output_directory is required without an orbit_solution directory");
		if (!prepared)
			keys(input, {"kind", "state", "path"});
		if (prepared)
			result.input = std::move(*prepared);
		else if (kind == "fixed_state" && !input.contains("path"))
			result.input = readEntry(input.at("state"), 0);
		else if (kind == "state_samples" && !input.contains("state"))
			result.input = EntryFile{resolve(base, input.at("path"))};
		else
			throw std::invalid_argument(
			    "StrewnCore inputs: fixed_state or state_samples; use StrewnOrbit for orbit_solution");
		const auto& environment = value.at("environment");
		const auto mode = environment.at("mode").get<std::string>();
		if (mode == "idealized") {
			keys(environment, {"mode", "ground_height_wgs84_m"});
			result.environment = IdealizedEnvironment{environment.value("ground_height_wgs84_m", 0.0)};
		} else if (mode == "file") {
			keys(environment, {"mode", "path"});
			result.environment = EnvironmentFile{resolve(base, environment.at("path"))};
		} else if (mode == "remote") {
			keys(environment, {"mode", "weather_url", "terrain_url", "timeout_seconds", "max_download_bytes"});
			const auto timeout = count(environment, "timeout_seconds", 60);
			const auto limit = count(environment, "max_download_bytes", 64 * 1024 * 1024);
			if (timeout == 0 || timeout > 3600 || limit == 0 || limit > 1024ULL * 1024 * 1024)
				throw std::invalid_argument("Invalid remote environment download limits");
			result.environment = EnvironmentUrls{environment.at("weather_url").get<std::string>(),
			                                     environment.at("terrain_url").get<std::string>(),
			                                     static_cast<unsigned>(timeout), static_cast<std::size_t>(limit)};
		} else if (mode == "automatic") {
			keys(environment, {"mode", "margin_km", "terrain_step_deg", "bounds_deg", "cycle_mjd_utc", "gfs_file",
			                   "timeout_seconds", "max_download_bytes", "offline", "on_unavailable"});
			AutomaticEnvironment options;
			options.marginKm = environment.value("margin_km", options.marginKm);
			options.terrainStepDeg = environment.value("terrain_step_deg", options.terrainStepDeg);
			if (environment.contains("bounds_deg"))
				options.boundsDeg = numbers<4>(environment.at("bounds_deg"));
			if (environment.contains("cycle_mjd_utc"))
				options.cycleMjdUtc = environment.at("cycle_mjd_utc").get<double>();
			if (environment.contains("gfs_file"))
				options.gfsFile = resolve(base, environment.at("gfs_file"));
			const auto timeout = count(environment, "timeout_seconds", options.timeoutSeconds);
			if (timeout > 3600)
				throw std::invalid_argument("Invalid automatic download timeout");
			options.timeoutSeconds = static_cast<unsigned>(timeout);
			options.maxDownloadBytes = count(environment, "max_download_bytes", options.maxDownloadBytes);
			options.offline = environment.value("offline", false);
			const auto fallback = environment.value("on_unavailable", "error");
			if (fallback != "error" && fallback != "idealized")
				throw std::invalid_argument("Unknown on_unavailable policy");
			options.idealizedFallback = fallback == "idealized";
			sources::validateAutomatic(options);
			result.environment = options;
		} else {
			throw std::invalid_argument("Supported environment modes: idealized, file, remote or automatic");
		}
		if (value.contains("physics")) {
			const auto& p = value.at("physics");
			keys(p, {"phase", "ram_pressure", "fragment_shape", "min_mass_kg", "luminous_speed_mps", "weibull_alpha",
			         "strength_cap_pa", "lateral_coefficient", "dust_at_every_breakup", "individual_fragments",
			         "bins_per_decade", "max_particles", "representatives_per_bin", "retain_mass_remainder"});
			auto& c = result.physics;
			const auto phase = p.value("phase", "entry");
			if (phase != "entry" && phase != "dark_flight")
				throw std::invalid_argument("Unknown flight phase");
			c.phase = phase == "entry" ? FlightPhase::Entry : FlightPhase::DarkFlight;
			const auto pressure = p.value("ram_pressure", "rho_v2");
			if (pressure != "rho_v2" && pressure != "half_rho_v2")
				throw std::invalid_argument("Unknown ram pressure convention");
			c.ramPressure = pressure == "rho_v2" ? RamPressure::RhoV2 : RamPressure::HalfRhoV2;
			const auto shape = p.value("fragment_shape", "sphere");
			if (shape == "sphere")
				c.fragmentShape = Shape::Sphere;
			else if (shape == "brick_area_only")
				c.fragmentShape = Shape::BrickAreaOnly;
			else if (shape == "cone_area_only")
				c.fragmentShape = Shape::ConeAreaOnly;
			else
				throw std::invalid_argument("Unknown shape; nonspherical models require explicit _area_only suffix");
			c.minMassKg = p.value("min_mass_kg", c.minMassKg);
			c.luminousSpeedMps = p.value("luminous_speed_mps", c.luminousSpeedMps);
			c.weibullAlpha = p.value("weibull_alpha", c.weibullAlpha);
			c.strengthCapPa = p.value("strength_cap_pa", c.strengthCapPa);
			c.lateralCoefficient = p.value("lateral_coefficient", c.lateralCoefficient);
			c.dustAtEveryBreakup = p.value("dust_at_every_breakup", c.dustAtEveryBreakup);
			c.retainMassRemainder = p.value("retain_mass_remainder", c.retainMassRemainder);
			c.maxParticles = count(p, "max_particles", c.maxParticles);
			c.representativesPerBin = count(p, "representatives_per_bin", c.representativesPerBin);
			if (p.contains("individual_fragments"))
				c.individualFragments = counts<std::size_t>(p.at("individual_fragments"));
			if (p.contains("bins_per_decade"))
				c.binsPerDecade = counts<unsigned>(p.at("bins_per_decade"));
		}
		if (value.contains("integration")) {
			const auto& p = value.at("integration");
			keys(p, {"relative_tolerance", "position_tolerance_m", "velocity_tolerance_mps", "mass_tolerance_kg",
			         "initial_step_s", "max_step_s", "min_step_s", "max_flight_time_s", "max_step_attempts"});
			auto& c = result.integration;
			c.relativeTolerance = p.value("relative_tolerance", c.relativeTolerance);
			c.positionToleranceM = p.value("position_tolerance_m", c.positionToleranceM);
			c.velocityToleranceMps = p.value("velocity_tolerance_mps", c.velocityToleranceMps);
			c.massToleranceKg = p.value("mass_tolerance_kg", c.massToleranceKg);
			c.initialStepS = p.value("initial_step_s", c.initialStepS);
			c.maxStepS = p.value("max_step_s", c.maxStepS);
			c.minStepS = p.value("min_step_s", c.minStepS);
			c.maxFlightTimeS = p.value("max_flight_time_s", c.maxFlightTimeS);
			c.maxStepAttempts = count(p, "max_step_attempts", c.maxStepAttempts);
		}
		if (value.contains("monte_carlo")) {
			const auto& p = value.at("monte_carlo");
			keys(p, {"seed", "runs", "workers", "parameters", "convergence", "run_timeout_s"});
			auto& c = result.monteCarlo;
			c.seed = count(p, "seed", c.seed);
			c.runs = count(p, "runs", c.runs);
			c.runTimeoutS = p.value("run_timeout_s", c.runTimeoutS);
			const auto workerCount = count(p, "workers", c.workers);
			if (workerCount > 1024)
				throw std::invalid_argument("workers must be <= 1024");
			c.workers = static_cast<unsigned>(workerCount);
			if (p.contains("convergence") && !p.at("convergence").is_null()) {
				const auto& v = p.at("convergence");
				keys(v, {"threshold", "batch_size", "consecutive_checks", "center_deg", "half_size_km", "cell_m",
				         "kernel_sigma_m"});
				ConvergenceConfig convergence;
				convergence.centerDeg = numbers<2>(v.at("center_deg"));
				convergence.threshold = v.value("threshold", convergence.threshold);
				convergence.batchSize = count(v, "batch_size", convergence.batchSize);
				convergence.consecutiveChecks = count(v, "consecutive_checks", convergence.consecutiveChecks);
				if (v.contains("half_size_km"))
					convergence.halfSizeKm = numbers<2>(v.at("half_size_km"));
				convergence.cellM = v.value("cell_m", convergence.cellM);
				convergence.kernelSigmaM = v.value("kernel_sigma_m", convergence.kernelSigmaM);
				detail::validateConvergence(convergence);
				c.convergence = convergence;
			}
			if (p.contains("parameters")) {
				const auto& v = p.at("parameters");
				keys(v, {"size", "density_kg_m3", "strong_strength_pa", "weak_strength_pa", "weak_probability",
				         "dust_fraction", "largest_fraction", "beta", "wind_scale", "drag_scale", "ablation_kg_j"});
				auto& r = c.parameters;
				r.densityKgM3 = range(v, "density_kg_m3", r.densityKgM3);
				r.strongStrengthPa = range(v, "strong_strength_pa", r.strongStrengthPa);
				r.weakStrengthPa = range(v, "weak_strength_pa", r.weakStrengthPa);
				r.weakProbability = v.value("weak_probability", r.weakProbability);
				r.dustFraction = range(v, "dust_fraction", r.dustFraction);
				r.largestFraction = range(v, "largest_fraction", r.largestFraction);
				r.beta = range(v, "beta", r.beta);
				r.windScale = range(v, "wind_scale", r.windScale);
				r.dragScale = range(v, "drag_scale", r.dragScale);
				r.ablationKgJ = range(v, "ablation_kg_j", r.ablationKgJ);
				if (v.contains("size")) {
					const auto& size = v.at("size");
					keys(size, {"mass_kg", "diameter_m", "absolute_magnitude", "albedo"});
					const int modes = static_cast<int>(size.contains("mass_kg")) +
					                  static_cast<int>(size.contains("diameter_m")) +
					                  static_cast<int>(size.contains("absolute_magnitude"));
					if (modes != 1 || (size.contains("albedo") && !size.contains("absolute_magnitude")))
						throw std::invalid_argument("Choose exactly one mass, diameter or H/albedo size model");
					if (size.contains("mass_kg"))
						r.size = MassSize{range(size, "mass_kg", {})};
					else if (size.contains("diameter_m"))
						r.size = DiameterSize{range(size, "diameter_m", {})};
					else {
						const auto albedo = numbers<3>(size.at("albedo"));
						r.size = PhotometricSize{size.at("absolute_magnitude").get<double>(), albedo[0], albedo[1],
						                         albedo[2]};
					}
				}
			}
		}
		validate(result.physics, result.integration);
		validate(result.monteCarlo);
		return result;
	}

	JobSpec loadJob(const std::filesystem::path& path) { return loadJobImpl(path, std::nullopt); }

	JobSpec loadJob(const std::filesystem::path& path, EntrySamples prepared)
	{
		return loadJobImpl(path, std::move(prepared));
	}

}  // namespace strewn
