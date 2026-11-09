#include "orbital_core/hazard/ManifoldOfVariations.h"
#include "orbital_core/TimeScale.h"
#include "fit/FitObservations.h"
#include "IntervalMcAssets.h"
#include "MovAssets.h"
#include <nlohmann/json.hpp>
#include <fstream>

namespace orbital_core::hazard {
	namespace {
		using Json = nlohmann::json;

		template <std::size_t N>
		std::string asset(const unsigned char (&v)[N])
		{
			return {reinterpret_cast<const char*>(v), N};
		}

		void replace(std::string& text, const std::string& key, const std::string& value)
		{
			std::size_t p = 0;
			while ((p = text.find(key, p)) != std::string::npos) {
				text.replace(p, key.size(), value);
				p += value.size();
			}
		}

		Json optional(std::optional<double> v) { return v ? Json(*v) : Json(nullptr); }

		Json vector(const Vec3D& v) { return {v.x, v.y, v.z}; }

		Json crossing(const std::optional<BoundaryCrossing>& entry)
		{
			if (!entry)
				return nullptr;
			return {{"epoch", entry->mjdTdb},
			        {"latitude_deg", entry->latitudeDeg},
			        {"longitude_deg", entry->longitudeDeg},
			        {"relative_position_km", vector(entry->relativePositionKm)},
			        {"relative_velocity_km_s", vector(entry->relativeVelocityKmS)}};
		}

		Json nodes(const std::vector<MovNode>& values)
		{
			auto rows = Json::array();
			for (const auto& n : values) {
				rows.push_back({{"row", n.row},
				                {"column", n.column},
				                {"status", static_cast<int>(n.status)},
				                {"attributable", n.attributable},
				                {"state", n.state},
				                {"chi_square", n.chiSquare},
				                {"chi", n.chi},
				                {"jacobian_mov", n.jacobianMov},
				                {"jacobian_sampling", n.jacobianSampling},
				                {"weight", n.weight},
				                {"H", optional(n.absoluteMagnitude)},
				                {"outer_boundary", n.outerBoundary},
				                {"entry_before_interval", n.entryBeforeInterval},
				                {"entry", crossing(n.entry)},
				                {"error", n.error}});
			}
			return rows;
		}
	}  // namespace

	void saveManifoldOfVariations(const MovResult& result, std::span<const ObsNode> observations,
	                              const FitConfig& settings, const std::filesystem::path& directory,
	                              const std::string& designation)
	{
		const auto& c = result.config;
		const auto& f = settings.nbody;
		Json data = {
		    {"schema", "orbital-mov-1"},
		    {"object", designation},
		    {"time_scale", "MJD TDB"},
		    {"state_frame", "heliocentric ecliptic J2000; AU, AU/day"},
		    {"attributable_units", "rad, rad, rad/day, rad/day, AU, AU/day"},
		    {"sampling", result.sampling == MovSampling::Cobweb ? "cobweb" : "grid"},
		    {"preliminary_log_range", result.preliminaryLogRange},
		    {"log_range", result.logRange},
		    {"density",
		     "exp(-chi^2/2) * sqrt(det(I + dA/drho^T dA/drho)) * sampling cell Jacobian; formal normal inverse"},
		    {"probability_scope",
		     "Conditional on sampled chi/H domain and resolved constrained fits; not a binomial confidence interval"},
		    {"epoch", result.epochMjdTdb},
		    {"attributable", result.attributable},
		    {"observer", result.observerCode},
		    {"curvature_snr", result.curvatureSnr},
		    {"start", c.startMjdTdb},
		    {"end", c.endMjdTdb},
		    {"target_name", longLabel(c.target)},
		    {"altitude_km", c.altitudeKm},
		    {"body_fixed_frame", result.bodyFixedFrame},
		    {"body_fixed_model", result.bodyFixedModel},
		    {"preliminary_nodes", c.preliminaryNodes},
		    {"range_nodes", c.rangeNodes},
		    {"rate_nodes", c.rateNodes},
		    {"sigma_max", c.sigmaMax},
		    {"H_max", c.absoluteMagnitudeMax},
		    {"workers", c.workers},
		    {"integrator", "RA15"},
		    {"forces",
		     {{"planet_mask", f.planetMask},
		      {"step_days", f.stepSize},
		      {"sb441", f.useSb441},
		      {"target_naif_id", f.targetNaifId},
		      {"interior_mass_policy", static_cast<int>(f.inactiveInteriorMassPolicy)}}},
		    {"weight_profile", static_cast<int>(settings.weightProfile)},
		    {"observer_sigmas", settings.obsWeighting.useObserverSigmas},
		    {"timing_sigma_default_seconds", settings.obsWeighting.timingSigmaDefaultSeconds},
		    {"cancelled", result.cancelled},
		    {"warnings", result.warnings},
		    {"outer_boundary_weight", result.outerBoundaryWeight},
		    {"probability_lower", optional(result.probabilityLower)},
		    {"probability_upper", optional(result.probabilityUpper)},
		    {"status_codes",
		     {"unstarted", "outside_region", "fit_failed", "rejected", "earth_bound", "complete", "propagation_failed",
		      "cancelled"}},
		    {"preliminary", nodes(result.preliminary)},
		    {"nodes", nodes(result.nodes)},
		    {"observations", Json::array()}};
		data["nominal"] = nullptr;
		if (result.nominal) {
			const auto& nominal = *result.nominal;
			data["nominal"] = {{"status", nominal.status == NominalPropagationStatus::Complete ? "complete"
			                             : nominal.status == NominalPropagationStatus::Cancelled ? "cancelled"
			                                                                        : "propagation_failed"},
			                   {"epoch", nominal.epochMjdTdb},
			                   {"state", nominal.state},
			                   {"entry_before_interval", nominal.entryBeforeInterval},
			                   {"entry", crossing(nominal.entry)},
			                   {"error", nominal.error}};
		}
		OrbitOptimizerConfig cfg;
		cfg.nbody = f;
		cfg.useWeights = settings.useWeights;
		cfg.weightProfile = settings.weightProfile;
		cfg.obsWeighting = settings.obsWeighting;
		const auto prepared = fit::prepareFitObservations(observations, cfg);
		for (std::size_t i = 0; i < observations.size(); ++i) {
			const auto& o = observations[i];
			const auto& p = prepared.at(i);
			data["observations"].push_back({{"epoch", p.mjdTdb},
			                                {"station", o.obsCode},
			                                {"ra_rad", o.RA},
			                                {"dec_rad", o.DEC},
			                                {"bias_ra_arcsec", o.biasRaArcsec},
			                                {"bias_dec_arcsec", o.biasDecArcsec},
			                                {"sigma_ra_rad", p.sigmaRa},
			                                {"sigma_dec_rad", p.sigmaDec},
			                                {"correlation", p.sigmaCorr},
			                                {"time_sigma_days", p.timeSigmaDays},
			                                {"earth_position_au", vector(p.earthR)},
			                                {"observer_offset_ecliptic_au", vector(p.observerOffsetEclAu)},
			                                {"mag", o.MAG},
			                                {"band", std::string(1, o.magBand)},
			                                {"mag_sigma", o.rmsMag},
			                                {"longitude", o.longitude},
			                                {"cos_phi", o.cosPhi},
			                                {"sin_phi", o.sinPhi}});
		}
		if (std::filesystem::exists(directory))
			throw std::runtime_error("MOV archive directory already exists");
		std::filesystem::create_directories(directory);
		const auto temporary = directory / "result.json.partial";
		std::ofstream out(temporary, std::ios::binary);
		out.exceptions(std::ios::badbit | std::ios::failbit);
		out << data.dump(2) << '\n';
		out.close();
		std::filesystem::rename(temporary, directory / "result.json");
	}

	void renderManifoldOfVariationsReport(const std::filesystem::path& directory)
	{
		std::ifstream input(directory / "result.json", std::ios::binary);
		if (!input)
			throw std::runtime_error("Cannot read MOV archive");
		Json data = Json::parse(input);
		if (data.at("schema") != "orbital-mov-1" || data.at("time_scale") != "MJD TDB")
			throw std::runtime_error("Unsupported MOV archive");
		for (const auto* key : {"epoch", "start", "end"})
			data[std::string(key) + "_utc"] = time::fromMjdTdb(data.at(key).get<double>(), time::Scale::UTC);
		for (auto& n : data.at("nodes"))
			if (n.at("entry").is_object()) {
				auto& e = n.at("entry");
				e["epoch_utc"] = time::fromMjdTdb(e.at("epoch").get<double>(), time::Scale::UTC);
			}
		std::string html = asset(mov_assets::html);
		replace(html, "__LEAFLET_CSS__", std::string(report_assets::leafletCss));
		replace(html, "__REPORT_CSS__", asset(report_assets::style));
		replace(html, "__LEAFLET_JS__", asset(report_assets::leaflet));
		replace(html, "__GRATICULE_JS__", asset(report_assets::graticule));
		replace(html, "__MAP_HELPERS__", asset(report_assets::script));
		replace(html, "__VIEWER_JS__", asset(mov_assets::script));
		replace(html, "__LICENSES__", std::string(report_assets::license));
		auto json = data.dump();
		replace(json, "<", "\\u003c");
		replace(html, "__DATA__", json);
		std::ofstream out(directory / "report.html", std::ios::binary);
		out.exceptions(std::ios::badbit | std::ios::failbit);
		out << html;
		out.close();
	}
}  // namespace orbital_core::hazard
