#include "JobFiles.h"
#include "VectorMath.h"

#include <algorithm>
#include <iomanip>
#include <locale>
#include <stdexcept>
#include <type_traits>

namespace strewn::jobfiles {

	std::string pathText(const std::filesystem::path& path)
	{
		const auto value = path.u8string();
		return {value.begin(), value.end()};
	}

	Json readJson(const std::filesystem::path& path)
	{
		std::ifstream stream(path, std::ios::binary);
		if (!stream)
			throw std::runtime_error("Cannot read " + pathText(path));
		return Json::parse(stream);
	}

	void writeJson(const std::filesystem::path& path, const Json& value)
	{
		std::ofstream stream(path, std::ios::binary | std::ios::trunc);
		stream.exceptions(std::ios::badbit | std::ios::failbit);
		stream << value.dump(2) << '\n';
	}

	void keys(const Json& object, std::initializer_list<std::string_view> allowed)
	{
		if (!object.is_object())
			throw std::invalid_argument("Expected a JSON object");
		for (const auto& [key, value] : object.items())
			if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
				throw std::invalid_argument("Unknown configuration key: " + key);
	}

	std::uint64_t count(const Json& object, const char* key, std::uint64_t fallback)
	{
		if (!object.contains(key))
			return fallback;
		const auto& value = object.at(key);
		if (!value.is_number_integer() || (!value.is_number_unsigned() && value.get<std::int64_t>() < 0))
			throw std::invalid_argument(std::string(key) + " must be a nonnegative integer");
		return value.get<std::uint64_t>();
	}

	Json jsonRange(Range value) { return Json::array({value.min, value.max}); }

	EntryState readEntry(const Json& value, std::uint64_t index, bool legacy)
	{
		EntryState result;
		if (legacy) {
			result.positionEcefM = numbers<3>(value.at("r_itrs"));
			result.velocityEcefMps = numbers<3>(value.at("v_itrs"));
			result.epochMjdUtc = value.at("jd_utc").get<double>() - 2400000.5;
		} else {
			keys(value, {"position_ecef_m", "velocity_ecef_mps", "epoch_mjd_utc", "epoch_jd_utc", "weight", "id",
			             "latitude_deg", "longitude_deg", "height_wgs84_m", "velocity_enu_mps"});
			if (value.contains("epoch_mjd_utc") == value.contains("epoch_jd_utc"))
				throw std::invalid_argument("Entry requires exactly one epoch_mjd_utc or legacy epoch_jd_utc");
			result.epochMjdUtc = value.contains("epoch_mjd_utc") ? value.at("epoch_mjd_utc").get<double>()
			                                                     : value.at("epoch_jd_utc").get<double>() - 2400000.5;
			const bool cartesian = value.contains("position_ecef_m") || value.contains("velocity_ecef_mps");
			const bool geodeticInput = value.contains("latitude_deg") || value.contains("longitude_deg") ||
			                           value.contains("height_wgs84_m") || value.contains("velocity_enu_mps");
			if (cartesian == geodeticInput)
				throw std::invalid_argument("Entry requires exactly one ECEF or geodetic/ENU representation");
			if (cartesian) {
				result.positionEcefM = numbers<3>(value.at("position_ecef_m"));
				result.velocityEcefMps = numbers<3>(value.at("velocity_ecef_mps"));
			} else {
				const GeodeticPosition position{value.at("latitude_deg").get<double>() * detail::deg,
				                                value.at("longitude_deg").get<double>() * detail::deg,
				                                value.at("height_wgs84_m").get<double>()};
				result.positionEcefM = ecefPosition(position);
				result.velocityEcefMps =
				    enuToEcef(numbers<3>(value.at("velocity_enu_mps")), position.latitudeRad, position.longitudeRad);
			}
		}
		result.weight = value.value("weight", 1.0);
		result.id = count(value, "id", index);
		return result;
	}

	Json entryJson(const EntryState& entry)
	{
		return {{"id", entry.id},
		        {"epoch_mjd_utc", entry.epochMjdUtc},
		        {"position_ecef_m", entry.positionEcefM},
		        {"velocity_ecef_mps", entry.velocityEcefMps},
		        {"weight", entry.weight}};
	}

	const char* shapeName(Shape shape)
	{
		switch (shape) {
			case Shape::Sphere:
				return "sphere";
			case Shape::BrickAreaOnly:
				return "brick_area_only";
			case Shape::ConeAreaOnly:
				return "cone_area_only";
		}
		throw std::invalid_argument("Unknown shape");
	}

	const char* runStatusName(RunStatus status)
	{
		switch (status) {
			case RunStatus::Complete:
				return "complete";
			case RunStatus::ParticleLimit:
				return "particle_limit";
			case RunStatus::Timeout:
				return "timeout";
			case RunStatus::Cancelled:
				return "cancelled";
			case RunStatus::Failed:
				return "failed";
			case RunStatus::WallTimeout:
				return "wall_timeout";
		}
		throw std::invalid_argument("Unknown run status");
	}

	Json parametersJson(const RunParameters& p)
	{
		return {{"mass_kg", p.massKg},
		        {"density_kg_m3", p.densityKgM3},
		        {"strength_pa", p.strengthPa},
		        {"dust_fraction", p.dustFraction},
		        {"largest_fraction", p.largestFraction},
		        {"beta", p.beta},
		        {"wind_scale", p.windScale},
		        {"drag_scale", p.dragScale},
		        {"ablation_kg_j", p.ablationKgJ}};
	}

	Json physicsJson(const PhysicsConfig& p)
	{
		return {{"phase", p.phase == FlightPhase::Entry ? "entry" : "dark_flight"},
		        {"ram_pressure", p.ramPressure == RamPressure::RhoV2 ? "rho_v2" : "half_rho_v2"},
		        {"fragment_shape", shapeName(p.fragmentShape)},
		        {"min_mass_kg", p.minMassKg},
		        {"luminous_speed_mps", p.luminousSpeedMps},
		        {"weibull_alpha", p.weibullAlpha},
		        {"strength_cap_pa", p.strengthCapPa},
		        {"lateral_coefficient", p.lateralCoefficient},
		        {"dust_at_every_breakup", p.dustAtEveryBreakup},
		        {"individual_fragments", p.individualFragments},
		        {"retain_mass_remainder", p.retainMassRemainder},
		        {"bins_per_decade", p.binsPerDecade},
		        {"representatives_per_bin", p.representativesPerBin},
		        {"max_particles", p.maxParticles}};
	}

	Json integrationJson(const IntegrationOptions& p)
	{
		return {{"relative_tolerance", p.relativeTolerance},
		        {"position_tolerance_m", p.positionToleranceM},
		        {"velocity_tolerance_mps", p.velocityToleranceMps},
		        {"mass_tolerance_kg", p.massToleranceKg},
		        {"initial_step_s", p.initialStepS},
		        {"max_step_s", p.maxStepS},
		        {"min_step_s", p.minStepS},
		        {"max_flight_time_s", p.maxFlightTimeS},
		        {"max_step_attempts", p.maxStepAttempts}};
	}

	Json samplingJson(const MonteCarloConfig& mc)
	{
		const auto& p = mc.parameters;
		Json size = std::visit(
		    [](const auto& size) -> Json {
			    using T = std::decay_t<decltype(size)>;
			    if constexpr (std::is_same_v<T, MassSize>)
				    return {{"mass_kg", jsonRange(size.kg)}};
			    else if constexpr (std::is_same_v<T, DiameterSize>)
				    return {{"diameter_m", jsonRange(size.m)}};
			    else
				    return {{"absolute_magnitude", size.absoluteMagnitude},
				            {"albedo", {size.albedoMin, size.albedoMode, size.albedoMax}}};
		    },
		    p.size);
		Json convergence = nullptr;
		if (mc.convergence) {
			const auto& c = *mc.convergence;
			convergence = {{"threshold", c.threshold},
			               {"batch_size", c.batchSize},
			               {"consecutive_checks", c.consecutiveChecks},
			               {"center_deg", c.centerDeg},
			               {"half_size_km", c.halfSizeKm},
			               {"cell_m", c.cellM},
			               {"kernel_sigma_m", c.kernelSigmaM}};
		}
		return {{"seed", mc.seed},
		        {"runs", mc.runs},
		        {"workers", mc.workers},
		        {"run_timeout_s", mc.runTimeoutS},
		        {"convergence", convergence},
		        {"parameters",
		         {{"size", size},
		          {"density_kg_m3", jsonRange(p.densityKgM3)},
		          {"strong_strength_pa", jsonRange(p.strongStrengthPa)},
		          {"weak_strength_pa", jsonRange(p.weakStrengthPa)},
		          {"weak_probability", p.weakProbability},
		          {"dust_fraction", jsonRange(p.dustFraction)},
		          {"largest_fraction", jsonRange(p.largestFraction)},
		          {"beta", jsonRange(p.beta)},
		          {"wind_scale", jsonRange(p.windScale)},
		          {"drag_scale", jsonRange(p.dragScale)},
		          {"ablation_kg_j", jsonRange(p.ablationKgJ)}}}};
	}

	Json summaryJson(const RunSummary& summary)
	{
		Json history = Json::array();
		for (const auto& check : summary.convergenceHistory)
			history.push_back({{"runs", check.runs},
			                   {"js_count", check.jsCount ? Json(*check.jsCount) : Json(nullptr)},
			                   {"js_mass", check.jsMass ? Json(*check.jsMass) : Json(nullptr)},
			                   {"stable_checks", check.stableChecks}});
		return {
		    {"requested", summary.requested},
		    {"completed", summary.completed},
		    {"truncated", summary.truncated},
		    {"failed", summary.failed},
		    {"wall_timeouts", summary.wallTimeouts},
		    {"runs_with_impacts", summary.runsWithImpacts},
		    {"tracked_impacts", summary.trackedImpacts},
		    {"workers", summary.workers},
		    {"cancelled", summary.cancelled},
		    {"elapsed_s", summary.elapsedS},
		    {"stop_reason", summary.cancelled   ? "cancelled"
		                    : summary.converged ? "js_converged"
		                                        : "run_limit"},
		    {"convergence_evaluated", std::any_of(summary.convergenceHistory.begin(), summary.convergenceHistory.end(),
		                                          [](const auto& c) { return c.jsCount.has_value(); })},
		    {"convergence_blocked", summary.convergenceBlocked},
		    {"convergence_out_of_grid_impacts", summary.convergenceOutOfGridImpacts},
		    {"convergence_history", history}};
	}

	Environment readEnvironment(const Json& value)
	{
		if (value.at("schema_version") != 1)
			throw std::invalid_argument("Unsupported environment schema");
		if (value.at("height_datum") != "WGS84")
			throw std::invalid_argument("Environment heights must be WGS84 ellipsoidal metres");
		const auto& a = value.at("atmosphere");
		if (a.contains("sample_columns") && a.at("sample_columns") != Json::array({"density_kg_m3", "temperature_k",
		                                                                           "east_mps", "north_mps", "up_mps"}))
			throw std::invalid_argument("Atmosphere sample columns do not match schema version 1");
		AtmosphereGrid grid;
		grid.latitudesDeg = a.at("latitudes_deg").get<std::vector<double>>();
		grid.longitudesDeg = a.at("longitudes_deg").get<std::vector<double>>();
		grid.firstHeightM = a.at("first_height_m").get<double>();
		grid.heightStepM = a.at("height_step_m").get<double>();
		grid.heightCount = count(a, "height_count", 0);
		for (const auto& sample : a.at("samples")) {
			const auto row = numbers<5>(sample);
			grid.samples.push_back({row[0], row[1], {row[2], row[3], row[4]}});
		}
		const auto& t = value.at("terrain");
		TerrainGrid terrain;
		terrain.latitudesDeg = t.at("latitudes_deg").get<std::vector<double>>();
		terrain.longitudesDeg = t.at("longitudes_deg").get<std::vector<double>>();
		terrain.heightsWgs84M = t.at("heights_wgs84_m").get<std::vector<double>>();
		return Environment(std::move(grid), std::move(terrain));
	}

	Json environmentJson(const Environment& environment)
	{
		const auto& a = environment.atmosphereGrid();
		const auto& t = environment.terrainGrid();
		Json rows = Json::array();
		for (const auto& sample : a.samples)
			rows.push_back({sample.densityKgM3, sample.temperatureK, sample.windEnuMps[0], sample.windEnuMps[1],
			                sample.windEnuMps[2]});
		return {{"schema_version", 1},
		        {"height_datum", "WGS84"},
		        {"atmosphere",
		         {{"latitudes_deg", a.latitudesDeg},
		          {"longitudes_deg", a.longitudesDeg},
		          {"first_height_m", a.firstHeightM},
		          {"height_step_m", a.heightStepM},
		          {"height_count", a.heightCount},
		          {"sample_columns", {"density_kg_m3", "temperature_k", "east_mps", "north_mps", "up_mps"}},
		          {"samples", rows}}},
		        {"terrain",
		         {{"latitudes_deg", t.latitudesDeg},
		          {"longitudes_deg", t.longitudesDeg},
		          {"heights_wgs84_m", t.heightsWgs84M}}}};
	}

	ResultArchive::ResultArchive(const std::filesystem::path& directory)
	    : m_impacts(directory / "impacts.csv"),
	      m_runs(directory / "runs.jsonl"),
	      m_breakups(directory / "breakups.jsonl")
	{
		for (auto* stream : {&m_impacts, &m_runs, &m_breakups})
			stream->exceptions(std::ios::badbit | std::ios::failbit);
		m_impacts.imbue(std::locale::classic());
		m_impacts
		    << std::setprecision(17)
		    << "run_id,entry_id,fragment_id,parent_id,latitude_deg,longitude_deg,height_wgs84_m,mass_kg,multiplicity,"
		       "speed_mps,time_since_entry_s,entry_mjd_utc,generation,run_weight,run_status,"
		       "x_ecef_m,y_ecef_m,z_ecef_m,vx_ecef_mps,vy_ecef_mps,vz_ecef_mps\n";
	}

	void ResultArchive::append(const RunResult& result)
	{
		for (const auto& impact : result.impacts) {
			const auto& f = impact.fragment;
			m_impacts << result.runId << ',' << result.entry.id << ',' << f.id << ',' << f.parentId << ','
			          << impact.latitudeDeg << ',' << impact.longitudeDeg << ',' << impact.heightWgs84M << ','
			          << f.state.massKg << ',' << f.multiplicity << ',' << detail::norm(f.state.velocityEcefMps) << ','
			          << f.timeS << ',' << result.entry.epochMjdUtc << ',' << f.generation << ",1,"
			          << runStatusName(result.status);
			for (double value : f.state.positionEcefM)
				m_impacts << ',' << value;
			for (double value : f.state.velocityEcefMps)
				m_impacts << ',' << value;
			m_impacts << '\n';
		}
		for (const auto& breakup : result.breakups) {
			const auto& p = breakup.parent;
			m_breakups << Json({{"run_id", result.runId},
			                    {"fragment_id", p.id},
			                    {"parent_id", p.parentId},
			                    {"time_s", p.timeS},
			                    {"generation", p.generation},
			                    {"mass_kg", p.state.massKg},
			                    {"multiplicity", p.multiplicity},
			                    {"height_wgs84_m", breakup.heightWgs84M},
			                    {"groups", breakup.groups},
			                    {"strength_pa", p.strengthPa},
			                    {"position_ecef_m", p.state.positionEcefM},
			                    {"velocity_ecef_mps", p.state.velocityEcefMps}})
			                  .dump()
			           << '\n';
		}
		const auto& b = result.mass;
		m_runs << Json({{"run_id", result.runId},
		                {"entry_id", result.entry.id},
		                {"run_weight", 1.0},
		                {"status", runStatusName(result.status)},
		                {"error", result.error},
		                {"particles", result.particles},
		                {"parameters", parametersJson(result.parameters)},
		                {"mass_balance_kg",
		                 {{"initial", b.initialKg},
		                  {"landed", b.landedKg},
		                  {"ablated", b.ablatedKg},
		                  {"dust", b.dustKg},
		                  {"residual", b.residualKg},
		                  {"below_cutoff", b.belowCutoffKg},
		                  {"unresolved", b.unresolvedKg}}},
		                {"integration",
		                 {{"evaluations", result.integration.evaluations},
		                  {"accepted_steps", result.integration.acceptedSteps},
		                  {"rejected_steps", result.integration.rejectedSteps}}}})
		              .dump()
		       << '\n';
		m_impacts.flush();
		m_breakups.flush();
		m_runs.flush();
	}


}  // namespace strewn::jobfiles
