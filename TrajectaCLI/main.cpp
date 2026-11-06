#include "orbital_core/ObservationSet.h"
#include <utility>
#include "hazard/LovCoordinateChart.h"
#include "NBodyPropagator.h"
#include "orbital_core/Constants.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Exporter.h"
#include "orbital_core/FitPipeline.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/CloseApproach.h"
#include "orbital_core/hazard/LineOfVariationImpact.h"
#include "orbital_core/forces/ForceAuditControl.h"
#include "orbital_core/hazard/ImpactRiskRmsObjective.h"
#include "orbital_core/hazard/ObsMonteCarlo.h"
#include "orbital_core/hazard/MonteCarloArchive.h"
#include "orbital_core/Photometry.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <nlohmann/json.hpp>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace {

	constexpr double kDefaultNBodyStepDays = 10.0;
	constexpr double kMinSaneNBodyStepDays = 1.0e-6;

	struct CliOptions {
		std::filesystem::path input;
		std::filesystem::path residualsPath;
		std::filesystem::path lovInputPath;
		orbital_core::FitConfig cfg;
		orbital_core::FitPipelineOptions pipeline;
		orbital_core::hazard::CloseApproachConfig caCfg;
		IntegratorKind surveyIntegrator = IntegratorKind::IAS15;
		bool debias = false;
		bool verbose = true;
		bool closeApproachScan = false;
		bool lovImpact = false;
		orbital_core::hazard::LineOfVariationImpactConfig lovCfg;
		int obsMcClones = 0;
		std::uint64_t obsMcSeed = 42;
		double obsMcAltitudeKm = 0.0;
		std::filesystem::path obsMcOutput;
		double outputEpochMjdOverride = 0.0;
		std::vector<double> stateAtMjd;
	};

	[[nodiscard]] std::string lower(std::string_view text)
	{
		std::string out(text);
		std::transform(out.begin(), out.end(), out.begin(),
		               [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
		return out;
	}

	[[nodiscard]] double parseDouble(const std::string& value, std::string_view name)
	{
		std::size_t pos = 0;
		const double parsed = std::stod(value, &pos);
		if (pos != value.size())
			throw std::invalid_argument("invalid value for " + std::string(name) + ": " + value);
		return parsed;
	}

	[[nodiscard]] int parseInt(const std::string& value, std::string_view name)
	{
		std::size_t pos = 0;
		const int parsed = std::stoi(value, &pos, 0);
		if (pos != value.size())
			throw std::invalid_argument("invalid value for " + std::string(name) + ": " + value);
		return parsed;
	}

	[[nodiscard]] uint16_t parseMask(const std::string& value)
	{
		const int parsed = parseInt(value, "--planet-mask");
		if (parsed < 0 || parsed > 0xffff)
			throw std::invalid_argument("--planet-mask must fit uint16");
		return static_cast<uint16_t>(parsed);
	}

	[[nodiscard]] double saneNBodyStepDays(double stepDays)
	{
		if (!std::isfinite(stepDays) || stepDays < kMinSaneNBodyStepDays)
			return kDefaultNBodyStepDays;
		return stepDays;
	}

	[[nodiscard]] orbital_core::ObservationSelectionPolicy parseObservationPolicy(const std::string& value)
	{
		const std::string v = lower(value);
		if (v == "extension" || v == "arc" || v == "arc-extension")
			return orbital_core::ObservationSelectionPolicy::ArcExtension;
		if (v == "all" || v == "all-original" || v == "global")
			return orbital_core::ObservationSelectionPolicy::AllOriginal;
		throw std::invalid_argument("unknown --obs-policy: " + value);
	}

	[[nodiscard]] orbital_core::AstrometricWeightProfile parseWeightProfile(const std::string& value)
	{
		const std::string v = lower(value);
		if (v == "vfcc17")
			return orbital_core::AstrometricWeightProfile::Vfcc17;
		if (v == "uniform")
			return orbital_core::AstrometricWeightProfile::Uniform;
		throw std::invalid_argument("unknown --weight-profile: " + value);
	}

	[[nodiscard]] orbital_core::NonGravModelKind parseNonGravModel(const std::string& value)
	{
		const std::string v = lower(value);
		if (v == "none")
			return orbital_core::NonGravModelKind::None;
		if (v == "h2o" || v == "marsden" || v == "marsdenh2o")
			return orbital_core::NonGravModelKind::MarsdenH2O;
		if (v == "co" || v == "marsdenco")
			return orbital_core::NonGravModelKind::MarsdenCO;
		if (v == "inverse-square" || v == "inversesquare" || v == "1/r2")
			return orbital_core::NonGravModelKind::InverseSquare;
		if (v == "yarkovsky")
			return orbital_core::NonGravModelKind::Yarkovsky;
		if (v == "dt" || v == "inverse-square-dt" || v == "inversesquaredt")
			return orbital_core::NonGravModelKind::InverseSquareDT;
		throw std::invalid_argument("unknown --non-grav-model: " + value);
	}

	[[nodiscard]] std::uint8_t parseNonGravSolveMask(const std::string& value)
	{
		std::string v = lower(value);
		v.erase(
		    std::remove_if(
		        v.begin(), v.end(),
		        [](char c) { return c == '-' || c == '_' || c == '/' || std::isspace(static_cast<unsigned char>(c)); }),
		    v.end());

		NgParameterMask mask = NgParameterMask::None;
		if (v == "default" || v == "model" || v == "none")
			mask = NgParameterMask::None;
		else if (v == "a2" || v == "yarkovsky")
			mask = NgParameterMask::A2;
		else if (v == "a12" || v == "a1a2")
			mask = NgParameterMask::A1 | NgParameterMask::A2;
		else if (v == "a123" || v == "a1a2a3")
			mask = NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3;
		else if (v == "a123dt" || v == "a1a2a3dt" || v == "all")
			mask = NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 | NgParameterMask::DT;
		else
			throw std::invalid_argument("unknown --ng-solve: " + value);

		return static_cast<std::uint8_t>(mask);
	}

	[[nodiscard]] orbital_core::hazard::Planet parseSurveyPlanet(const std::string& value)
	{
		const std::string v = lower(value);
		if (v == "mercury" || v == "merc")
			return orbital_core::hazard::Planet::Mercury;
		if (v == "venus")
			return orbital_core::hazard::Planet::Venus;
		if (v == "earth")
			return orbital_core::hazard::Planet::Earth;
		if (v == "mars")
			return orbital_core::hazard::Planet::Mars;
		if (v == "jupiter" || v == "juptr")
			return orbital_core::hazard::Planet::Jupiter;
		if (v == "saturn" || v == "satrn")
			return orbital_core::hazard::Planet::Saturn;
		if (v == "uranus" || v == "urnus")
			return orbital_core::hazard::Planet::Uranus;
		if (v == "neptune" || v == "neptn")
			return orbital_core::hazard::Planet::Neptune;
		if (v == "moon")
			return orbital_core::hazard::Planet::Moon;
		throw std::invalid_argument("unknown --survey-planet: " + value);
	}

	void printUsage(std::ostream& os)
	{
		os << "Usage: TrajectaCLI <mpc-file> [options]\n"
		   << "\n"
		   << "Options:\n"
		   << "  --obs-policy extension|all       Initial final-fit mask (default: all)\n"
		   << "  --weight-profile vfcc17|uniform Astrometric sigmas (uniform: 0.5 arcsec; default: vfcc17)\n"
		   << "  --epoch <mjd>                     Output epoch (TDB) for elements and sigmas\n"
		   << "  --time-sigma <sec>                Timing-uncertainty term: ADES rmsTime\n"
		   << "                                    where reported, this value elsewhere\n"
		   << "                                    (0 = reported only, 1 = customary\n"
		   << "                                    default; omit the flag to disable it)\n"
		   << "                                    (default: derived from the arc)\n"
		   << "  --observer-sigmas / --no-observer-sigmas\n"
		   << "                                    Weight by the observer-supplied ADES\n"
		   << "                                    rmsRA/rmsDec where present, table\n"
		   << "                                    elsewhere (default: off)\n"
		   << "  --debias / --no-debias            Apply catalog debias (default: on)\n"
		   << "  --weights / --no-weights          Use MPC weighting (default: on)\n"
		   << "  --sigma-clip / --no-sigma-clip    FitSession hysteretic sigma clip (default: on)\n"
		   << "  --sigma <N>                       Sigma clip threshold (default: 3)\n"
		   << "  --nbody / --no-nbody              N-body final fit (default: on)\n"
		   << "  --step <days>                     Initial integration step (default: 10)\n"
		   << "  --planet-mask <int>               Perturber mask (default: planets+Moon)\n"
		   << "  --sb441                           Enable SB441-N16 if available\n"
		   << "  --non-grav                        Enable non-grav fit\n"
		   << "  --non-grav-model <name>           h2o|co|inverse-square|yarkovsky|dt\n"
		   << "  --ng-solve <set>                  default|a12|a123|a123dt|a2\n"
		   << "  --init-a1 <AU/d^2> --init-a2 <AU/d^2> --init-a3 <AU/d^2>\n"
		   << "  --residuals <file>                Save final residual table\n"
		   << "  --dump-lov-input <file>           Full-precision nominal state and COV\n"
		   << "  --ca-scan                         Print nominal close approaches\n"
		   << "  --survey-horizon <years>          Close-approach horizon (default: 100)\n"
		   << "  --survey-integrator <name>        ra15|ias15 for CA/LOV/MC; fit unchanged (default: ias15)\n"
		   << "  --survey-lookback <years>         Lookback before fit epoch (default: 0)\n"
		   << "  --survey-planet <name>            Target body for CA scan (default: earth)\n"
		   << "  --survey-threshold <AU>           Report/filter radius/dmeacontr (default: 0.2)\n"
		   << "  --ca-encounter-radius <AU>        TP recording radius/dmea (default: 0.1)\n"
		   << "  --survey-step <days>              CA propagation sample step (default: 1)\n"
		   << "  --lov-impact                      Run LOV impact-monitor trace, print risk rows\n"
		   << "  --lov-scaling <on|off>            LOV coordinate scaling (default: on; independent of rescov)\n"
		   << "  --lov-sampling <mode>            fixed|probability (default: fixed)\n"
		   << "  --lov-prob-step <x>               Step for probability sampling (default: 1e-7)\n"
		   << "  --lov-sigma-step <x>              Fixed/max LOV sigma step (default: 0.01)\n"
		   << "  --obs-mc <N>                      State-covariance Monte Carlo with N clones\n"
		   << "  --obs-mc-seed <n>                 RNG seed for --obs-mc (default: 42)\n"
		   << "  --mc-altitude-km <h>              MC boundary height; 0 = legacy sphere\n"
		   << "  --mc-output <path>                Save Earth MC crossing states as JSON\n"
		   << "  --state-at <mjd>[,<mjd>...]       Propagate the fitted nominal to these\n"
		   << "                                    epochs (MJD TDB) and print the ECLM\n"
		   << "                                    J2000 state at full precision\n"
		   << "  --quiet                           Suppress progress logs\n"
		   << "  -h, --help                        Show this help\n";
	}

	[[nodiscard]] CliOptions parseArgs(int argc, char* argv[])
	{
		CliOptions opt;
		opt.cfg.weightProfile = orbital_core::AstrometricWeightProfile::Vfcc17;
		opt.debias = true;
		opt.cfg.nbody.enabled = true;
		opt.cfg.nbody.stepSize = kDefaultNBodyStepDays;
		opt.cfg.nbody.planetMask = static_cast<uint16_t>(PLANET_BODY_MASK | MOON_BODY_MASK);
		opt.cfg.nbody.useSb441 = false;
		opt.cfg.useWeights = true;
		opt.cfg.sigmaClip.enabled = true;
		opt.cfg.nonGrav.model = orbital_core::NonGravModelKind::MarsdenH2O;

		auto requireValue = [&](int& i, std::string_view name) -> std::string {
			if (i + 1 >= argc)
				throw std::invalid_argument(std::string(name) + " requires a value");
			return argv[++i];
		};

		for (int i = 1; i < argc; ++i) {
			const std::string arg = argv[i];
			if (arg == "-h" || arg == "--help") {
				printUsage(std::cout);
				std::exit(0);
			}
			if (!arg.empty() && arg[0] == '-') {
				if (arg == "--obs-policy") {
					opt.pipeline.observationSelectionPolicy = parseObservationPolicy(requireValue(i, arg));
				} else if (arg == "--weight-profile")
					opt.cfg.weightProfile = parseWeightProfile(requireValue(i, arg));
				else if (arg == "--epoch")
					opt.outputEpochMjdOverride = std::stod(requireValue(i, arg));
				else if (arg == "--observer-sigmas")
					opt.cfg.obsWeighting.useObserverSigmas = true;
				else if (arg == "--no-observer-sigmas")
					opt.cfg.obsWeighting.useObserverSigmas = false;
				else if (arg == "--time-sigma")
					opt.cfg.obsWeighting.timingSigmaDefaultSeconds = std::stod(requireValue(i, arg));
				else if (arg == "--debias")
					opt.debias = true;
				else if (arg == "--no-debias")
					opt.debias = false;
				else if (arg == "--weights")
					opt.cfg.useWeights = true;
				else if (arg == "--no-weights")
					opt.cfg.useWeights = false;
				else if (arg == "--sigma-clip")
					opt.cfg.sigmaClip.enabled = true;
				else if (arg == "--no-sigma-clip")
					opt.cfg.sigmaClip.enabled = false;
				else if (arg == "--sigma")
					opt.cfg.sigmaClip.thresholdNSigma = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--nbody")
					opt.cfg.nbody.enabled = true;
				else if (arg == "--no-nbody")
					opt.cfg.nbody.enabled = false;
				else if (arg == "--step")
					opt.cfg.nbody.stepSize = saneNBodyStepDays(parseDouble(requireValue(i, arg), arg));
				else if (arg == "--planet-mask")
					opt.cfg.nbody.planetMask = parseMask(requireValue(i, arg));
				else if (arg == "--sb441")
					opt.cfg.nbody.useSb441 = true;
				else if (arg == "--non-grav")
					opt.cfg.nonGrav.enabled = true;
				else if (arg == "--non-grav-model")
					opt.cfg.nonGrav.model = parseNonGravModel(requireValue(i, arg));
				else if (arg == "--ng-solve")
					opt.cfg.nonGrav.solveMask = parseNonGravSolveMask(requireValue(i, arg));
				else if (arg == "--init-a1") {
					opt.cfg.nonGrav.initialA1 = parseDouble(requireValue(i, arg), arg);
					opt.cfg.nonGrav.useInitialValues = true;
				} else if (arg == "--init-a2") {
					opt.cfg.nonGrav.initialA2 = parseDouble(requireValue(i, arg), arg);
					opt.cfg.nonGrav.useInitialValues = true;
				} else if (arg == "--init-a3") {
					opt.cfg.nonGrav.initialA3 = parseDouble(requireValue(i, arg), arg);
					opt.cfg.nonGrav.useInitialValues = true;
				} else if (arg == "--init-dt") {
					opt.cfg.nonGrav.initialDT = parseDouble(requireValue(i, arg), arg);
					opt.cfg.nonGrav.useInitialValues = true;
				} else if (arg == "--residuals")
					opt.residualsPath = requireValue(i, arg);
				else if (arg == "--dump-lov-input")
					opt.lovInputPath = requireValue(i, arg);
				else if (arg == "--ca-scan")
					opt.closeApproachScan = true;
				else if (arg == "--survey-horizon")
					opt.caCfg.horizonYears = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--survey-integrator") {
					const auto name = lower(requireValue(i, arg));
					if (name == "ra15")
						opt.surveyIntegrator = IntegratorKind::RA15;
					else if (name == "ias15")
						opt.surveyIntegrator = IntegratorKind::IAS15;
					else
						throw std::invalid_argument("unknown --survey-integrator: " + name);
				} else if (arg == "--survey-lookback")
					opt.caCfg.lookbackYears = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--survey-planet")
					opt.caCfg.planet = parseSurveyPlanet(requireValue(i, arg));
				else if (arg == "--survey-threshold")
					opt.caCfg.thresholdAu = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--ca-encounter-radius")
					opt.caCfg.encounterRadiusAu = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--survey-step")
					opt.caCfg.sampleStepDay = saneNBodyStepDays(parseDouble(requireValue(i, arg), arg));
				else if (arg == "--lov-impact")
					opt.lovImpact = true;
				else if (arg == "--lov-scaling") {
					const auto value = requireValue(i, arg);
					if (value != "on" && value != "off")
						throw std::runtime_error("--lov-scaling must be on or off");
					opt.lovCfg.coordinateScaling = value == "on";
				} else if (arg == "--lov-sampling") {
					const auto value = requireValue(i, arg);
					if (value != "fixed" && value != "probability")
						throw std::runtime_error("--lov-sampling must be fixed or probability");
					opt.lovCfg.probabilitySampling = value == "probability";
				} else if (arg == "--lov-prob-step")
					opt.lovCfg.probabilityStep = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--lov-sigma-step")
					opt.lovCfg.sigmaStepMax = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--obs-mc")
					opt.obsMcClones = static_cast<int>(parseDouble(requireValue(i, arg), arg));
				else if (arg == "--state-at") {
					const std::string list = requireValue(i, arg);
					std::size_t pos = 0;
					while (pos <= list.size()) {
						const std::size_t comma = list.find(',', pos);
						const std::string tok =
						    list.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
						if (!tok.empty())
							opt.stateAtMjd.push_back(parseDouble(tok, "--state-at"));
						if (comma == std::string::npos)
							break;
						pos = comma + 1;
					}
					if (opt.stateAtMjd.empty())
						throw std::invalid_argument("--state-at: no epochs parsed");
					std::sort(opt.stateAtMjd.begin(), opt.stateAtMjd.end());
				} else if (arg == "--obs-mc-seed")
					opt.obsMcSeed = static_cast<std::uint64_t>(std::stoull(requireValue(i, arg)));
				else if (arg == "--mc-altitude-km")
					opt.obsMcAltitudeKm = parseDouble(requireValue(i, arg), arg);
				else if (arg == "--mc-output")
					opt.obsMcOutput = requireValue(i, arg);
				else if (arg == "--quiet")
					opt.verbose = false;
				else if (arg.starts_with("-qmljsdebugger=")) {
				} else
					throw std::invalid_argument("unknown option: " + arg);
				continue;
			}

			if (!opt.input.empty())
				throw std::invalid_argument("multiple input files provided");
			opt.input = arg;
		}

		if (opt.input.empty())
			throw std::invalid_argument("missing input file");
		if (!std::isfinite(opt.obsMcAltitudeKm) || opt.obsMcAltitudeKm < 0.0)
			throw std::invalid_argument("--mc-altitude-km must be finite and nonnegative");
		if ((opt.obsMcAltitudeKm > 0.0 || !opt.obsMcOutput.empty()) && opt.obsMcClones <= 0)
			throw std::invalid_argument("MC altitude/output requires --obs-mc <N>");
		if (!opt.obsMcOutput.empty() && std::filesystem::exists(opt.obsMcOutput))
			throw std::invalid_argument("MC output already exists");
		return opt;
	}

	[[nodiscard]] orbital_core::FitPipelineResult runPipeline(const CliOptions& opt)
	{
		if (opt.verbose)
			std::cerr << "Loading input: " << opt.input.string() << '\n';
		auto observations = orbital_core::loadObservations(opt.input);
		auto& allObs = observations.all;
		if (opt.verbose) {
			std::cout << "Loaded raw observations: " << allObs.size() << " (IOD pool " << observations.iod.size()
			          << ")\n";
		}
		orbital_core::prepareObservations(observations, {.debias = opt.debias, .weightProfile = opt.cfg.weightProfile});
		if (opt.debias && opt.verbose) {
			const auto& preparation = *observations.preparation;
			std::cout << "AstrometryBias: debiased " << preparation.correctedAll << "/" << allObs.size()
			          << " observations; " << preparation.correctedIod << "/" << observations.iod.size()
			          << " IOD picks\n";
		}

		orbital_core::FitPipelineInput input{observations, std::nullopt};
		if (opt.outputEpochMjdOverride > 0.0) {
			input.outputEpochJdTdb = opt.outputEpochMjdOverride + 2400000.5;
			std::cout << "[Astrometry] output epoch forced to MJD " << std::fixed << std::setprecision(6)
			          << opt.outputEpochMjdOverride << " (JD " << *input.outputEpochJdTdb << ")\n"
			          << std::defaultfloat;
		}
		orbital_core::FitRunControl control;
		if (opt.verbose) {
			control.progress = [](const std::string& message) { std::cout << message << '\n'; };
		}
		auto result = orbital_core::FitPipeline{opt.cfg, opt.pipeline}.run(input, control);
		if (result.fit.cancelled)
			throw std::runtime_error("Cancelled");
		return result;
	}

	void printSummary(const orbital_core::FitPipelineResult& result)
	{
		std::cout << std::fixed << std::setprecision(12);
		std::cout << "IOD winner: " << result.initial.method << '\n';
		std::cout << "Loaded observations: " << result.observations.all.size() << '\n';
		if (result.extension && result.extension->observationCount > 0) {
			const auto& extension = *result.extension;
			std::cout << "ATT-EXT warm-start obs: " << extension.observationCount << '\n';
			std::cout << "ATT-EXT perturbers: configured=0x" << std::hex << extension.configuredMask << " effective=0x"
			          << extension.effectiveMask << " auto=0x" << extension.automaticMask << " blocked=0x"
			          << extension.blockedMask << std::dec << '\n';
		}
		std::cout << "Fit input observations: " << result.fitInputObservations << '\n';
		std::cout << "Final observations: " << result.fittedObservations().size() << '\n';
		std::cout << "RMS: " << std::setprecision(6) << result.fit.residuals.rms() << "\" ("
		          << result.fit.residuals.count() << " residual rows)\n";
		std::cout << "Cost: " << std::setprecision(12) << result.fit.cost << '\n';
		std::cout << "Iterations: " << result.fit.iterations
		          << " converged=" << (result.fit.converged ? "true" : "false") << '\n';
		if (result.fit.sigmaClipPassesRun > 0) {
			std::cout << "Sigma-clip passes: " << result.fit.sigmaClipPassesRun
			          << " rejected=" << result.fit.sigmaClipRejected << '\n';
		}

		const auto& orbit = result.fit.orbit;
		if (!orbit) {
			std::cout << "Final orbit: <null>\n";
			return;
		}

		std::cout << std::setprecision(15);
		std::cout << "a = " << orbit->getSemiMajorAxis() << " AU\n";
		std::cout << "e = " << orbit->getExcentricity() << '\n';
		std::cout << "i = " << Convertors::radiansToDegrees(orbit->getInclination()) << " deg\n";
		std::cout << "Omega = " << Convertors::radiansToDegrees(orbit->getLongitudeOfAscendingNode()) << " deg\n";
		std::cout << "omega = " << Convertors::radiansToDegrees(orbit->getArgumentOfPericentre()) << " deg\n";
		std::cout << "M = " << Convertors::radiansToDegrees(orbit->getMeanAnomaly()) << " deg\n";

		const NonGravModel& ng = result.fit.finalNonGrav;
		if (ng.enabled()) {
			std::cout << std::scientific << std::setprecision(6);
			std::cout << "A1 = " << ng.A1() << " AU/d^2\n";
			std::cout << "A2 = " << ng.A2() << " AU/d^2\n";
			std::cout << "A3 = " << ng.A3() << " AU/d^2\n";
			std::cout << std::fixed << std::setprecision(6);
			std::cout << "DT = " << ng.DT() << " d\n";
		}
	}

	[[nodiscard]] NBodyPropagator makeSurveyPropagator(const orbital_core::FitConfig& cfg,
	                                                   const orbital_core::FitResult& fit,
	                                                   IntegratorKind integrator = IntegratorKind::RA15,
	                                                   double encounterRadiusAu = Constants::EARTH_TARGET_PLANE_AU)
	{
		NBodyPropagator prop(integrator);
		prop.setStepSize(cfg.nbody.stepSize);
		prop.setPlanetMask(cfg.nbody.planetMask);
		prop.setUseSb441(cfg.nbody.useSb441);
		prop.setTargetNaifId(cfg.nbody.targetNaifId);
		prop.setInactiveInteriorMassPolicy(cfg.nbody.inactiveInteriorMassPolicy);
		prop.nonGrav() = fit.finalNonGrav;
		if (IEphemerisSource* eph = getDefaultEphemerisSource())
			prop = prop.forEncounterMonitoring(*eph, Constants::ENCOUNTER_STEP_POINTS, encounterRadiusAu);
		return prop;
	}

	void printCloseApproaches(const std::vector<orbital_core::hazard::CloseApproachEvent>& events)
	{
		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
		std::cout << std::fixed << std::setprecision(6);
		std::cout << "Nominal close approaches: " << events.size() << '\n';
		for (std::size_t i = 0; i < events.size(); ++i) {
			const auto& e = events[i];
			const double distKm = e.distAu * kAuKm;
			const double distLd = distKm / Constants::LUNAR_DISTANCE_KM;
			std::cout << "  ca[" << i << "]"
			          << " JD_TDB=" << std::fixed << std::setprecision(12) << e.jd << " UTC="
			          << Convertors::convertJDToGregorian(orbital_core::time::jdFromMjd(orbital_core::time::fromMjdTdb(
			                 orbital_core::time::mjdFromJd(e.jd), orbital_core::time::Scale::UTC)))
			          << " distAu=" << std::scientific << std::setprecision(15) << e.distAu << " distKm=" << std::fixed
			          << std::setprecision(6) << distKm << " distLD=" << std::setprecision(6) << distLd
			          << " vInfKmS=" << std::setprecision(3) << e.vInfKmS << '\n';
		}
	}

	void dumpLovInput(const std::filesystem::path& path, const orbital_core::FitResult& fit)
	{
		const Orbit_p orbit = fit.fitEpochOrbit ? fit.fitEpochOrbit : fit.orbit;
		const auto& covariance = fit.uncertainty.fitParameterCovariance;
		if (!orbit || covariance.rows() != 6 || covariance.cols() != 6 || fit.finalNonGrav.numFitParams() != 0 ||
		    !covariance.allFinite())
			throw std::runtime_error("--dump-lov-input requires a finite six-parameter fit");
		const bool cometary = orbital_core::cometary::usesCometaryElements(orbit->getExcentricity());
		const auto chart = cometary
		                       ? orbital_core::hazard::makeCometaryLovChart(*orbit, covariance, fit.finalNonGrav)
		                       : orbital_core::hazard::makeEquinoctialLovChart(*orbit, covariance, fit.finalNonGrav);
		if (!chart.valid)
			throw std::runtime_error(cometary ? "--dump-lov-input: invalid cometary chart"
			                                  : "--dump-lov-input: invalid equinoctial chart");
		if (std::filesystem::exists(path))
			throw std::runtime_error("--dump-lov-input refuses to overwrite " + path.string());
		std::ofstream out(path);
		if (!out)
			throw std::runtime_error("cannot open --dump-lov-input file");
		out.imbue(std::locale::classic());
		out << std::scientific << std::setprecision(17)
		    << (cometary ? "LOV_INPUT 2\nCOORDINATE_SYSTEM COT\nEPOCH_JD " : "LOV_INPUT 1\nEPOCH_JD ")
		    << orbit->getEpoch() << "\nCOVARIANCE_OBJECTIVE_SCALE " << fit.uncertainty.covarianceObjectiveScale;
		const Vec3D r = orbit->getR();
		const Vec3D v = orbit->getV() * Constants::GAUSS_K;
		out << "\nSTATE " << r.x << ' ' << r.y << ' ' << r.z << ' ' << v.x << ' ' << v.y << ' ' << v.z
		    << (cometary ? "\nCOT_RAD" : "\nEQU_RAD");
		for (int i = 0; i < 6; ++i)
			out << ' ' << chart.nominalCoordinates(i);
		out << '\n';
		for (int i = 0; i < 6; ++i) {
			out << "COV_PHYSICAL " << i;
			for (int j = 0; j < 6; ++j)
				out << ' ' << covariance(i, j);
			out << (cometary ? "\nCOV_COT_RAD " : "\nCOV_EQU_RAD ") << i;
			for (int j = 0; j < 6; ++j)
				out << ' ' << chart.coordinateCovariance(i, j);
			out << '\n';
		}
		out.flush();
		if (!out)
			throw std::runtime_error("failed writing --dump-lov-input file");
	}

	void printLovImpact(const Orbit& orbit, const orbital_core::FitResult& fit,
	                    std::span<const ObsNode> fitObservations, const orbital_core::FitConfig& fitConfig,
	                    const NBodyPropagator& propagator, orbital_core::hazard::LineOfVariationImpactConfig cfg)
	{
		namespace hz = orbital_core::hazard;
		const auto& unc = fit.uncertainty;
		const bool hasParamCov = unc.fitParameterCovariance.size() > 0;
		const bool hasStateCov = !unc.cartesianStateCovariance.isZero(0.0);
		if (!hasParamCov && !hasStateCov) {
			std::cout << "LOV impact: fit produced no covariance\n";
			return;
		}

		cfg.closeApproach.sampleStepDay = 1.0;
		if (!fitObservations.empty()) {
			cfg.impactRefinement.rmsScaling =
			    hz::makeImpactRiskRmsScaling(std::make_shared<Orbit>(orbit), fitObservations, fit.keptObsIndices,
			                                 fitConfig, fit.finalNonGrav, unc.covarianceObjectiveScale);
		}

		orbital_core::forces::audit::arm("lov");

		const std::vector<hz::ImpactRiskRow> rows =
		    hasParamCov
		        ? hz::traceLineOfVariationImpacts(orbit, unc.fitParameterCovariance, fit.finalNonGrav, propagator, cfg)
		        : hz::traceLineOfVariationImpacts(orbit, unc.cartesianStateCovariance, propagator, cfg);

		std::cout << "LOV impact rows: " << rows.size() << (hasParamCov ? " (parameter-cov)\n" : " (state-cov)\n");
		for (std::size_t i = 0; i < rows.size(); ++i) {
			const auto& res = rows[i].result;
			std::cout << "  row[" << i << "]"
			          << " jd=" << std::fixed << std::setprecision(3) << rows[i].cluster.meanJd
			          << " p=" << std::scientific << std::setprecision(6) << res.probabilityImpact
			          << " sigmaLov=" << std::fixed << std::setprecision(6) << res.sigmaLov << " gamma=" << res.gamma
			          << " vInf=" << std::setprecision(3) << res.vInfKmS << " minMissKm=" << std::setprecision(1)
			          << res.minMissKm << " cluster=" << rows[i].clusterIndex << " path=" << rows[i].pathIndex << "/"
			          << rows[i].pathCount << '\n';
		}
	}

	void printObsMc(const Orbit& orbit, const orbital_core::FitResult& fit, const NBodyPropagator& propagator,
	                const std::vector<orbital_core::hazard::CloseApproachEvent>& events, int numClones,
	                std::uint64_t seed, double altitudeKm, const std::filesystem::path& output,
	                const std::filesystem::path& solutionDirectory, const std::string& objectName)
	{
		namespace hz = orbital_core::hazard;
		if (events.empty()) {
			std::cout << "cov-MC: no close-approach events\n";
			return;
		}
		if (fit.uncertainty.cartesianStateCovariance.isZero(0.0)) {
			std::cout << "cov-MC: fit produced no state covariance\n";
			return;
		}

		hz::ObsMcConfig cfg;
		cfg.numClones = numClones;
		cfg.rngSeed = seed;
		cfg.altitudeKm = altitudeKm;

		const auto& fullCov = fit.uncertainty.fitParameterCovariance;
		const bool full = fullCov.rows() == 6 + fit.finalNonGrav.numFitParams() && fullCov.rows() == fullCov.cols();
		const hz::ObsMcResult mc =
		    full ? hz::runParameterCovarianceMonteCarlo(orbit, fullCov, fit.finalNonGrav, propagator, events, cfg)
		         : hz::runStateCovarianceMonteCarlo(orbit, fit.uncertainty.cartesianStateCovariance, propagator, events,
		                                            cfg);
		hz::MonteCarloArchiveInput archive;
		archive.solution.nominal = orbit;
		archive.solution.objectName = objectName;
		archive.solution.nonGrav = fit.finalNonGrav;
		archive.solution.absoluteMagnitude =
		    orbital_core::photometry::fit(orbital_core::photometry::observations(fit).data()).absoluteMagnitude;
		if (full)
			archive.solution.covariance = fullCov;
		else
			archive.solution.covariance = fit.uncertainty.cartesianStateCovariance;
		archive.kind = full ? hz::CovarianceSampling::Parameters : hz::CovarianceSampling::State;
		archive.sampling = cfg;
		archive.covarianceObjectiveScale = fit.uncertainty.covarianceObjectiveScale;
		const auto saved = hz::saveMonteCarloArchive(solutionDirectory, archive, propagator, events, mc);
		std::cout << "MC archive written to: " << saved.directory.string() << '\n';

		if (altitudeKm > 0.0)
			std::cout << "cov-MC boundary=" << altitudeKm
			          << " km above reference ellipsoid; counts/P refer to entry, not ground impact\n";

		if (!output.empty()) {
			using nlohmann::json;
			json data = {{"schema_version", 1},
			             {"frame", "ECEF"},
			             {"frame_model", "OrbitalCore IAU1976/1980 UT1/GAST; no polar motion"},
			             {"altitude_km", altitudeKm},
			             {"boundary", altitudeKm > 0.0 ? "reference_ellipsoid_altitude" : "planet_sphere"},
			             {"seed", seed},
			             {"requested", numClones},
			             {"completed", mc.numConverged},
			             {"failed", mc.numFailedFit},
			             {"cancelled", mc.cancelled},
			             {"conditioning", "Earth boundary crossings within selected encounter windows"},
			             {"entries", json::array()},
			             {"crossings", json::array()}};
			for (const auto& point : mc.corridor) {
				if (point.planet != hz::Planet::Earth)
					continue;
				const auto id = static_cast<std::uint64_t>(point.cloneIndex) * events.size() +
				                static_cast<std::uint64_t>(point.eventIndex);
				const auto r = point.rBodyFixedKm * 1000.0;
				const auto v = point.vBodyFixedKmS * 1000.0;
				data["entries"].push_back({{"id", id},
				                           {"epoch_mjd_utc", point.epochMjdUtc},
				                           {"position_ecef_m", {r.x, r.y, r.z}},
				                           {"velocity_ecef_mps", {v.x, v.y, v.z}},
				                           {"weight", 1.0}});
				data["crossings"].push_back({{"id", id},
				                             {"clone", point.cloneIndex},
				                             {"event", point.eventIndex},
				                             {"epoch_mjd_tdb", point.epochMjdTdb},
				                             {"latitude_deg", point.latitudeDeg},
				                             {"longitude_deg", point.longitudeDeg}});
			}
			std::ofstream stream(output);
			stream.exceptions(std::ios::failbit | std::ios::badbit);
			stream << data.dump(2) << '\n';
		}

		std::cout << "cov-MC clones=" << mc.numAttempted << " converged=" << mc.numConverged
		          << " anyImpact=" << mc.numAnyImpactClones << " aggregateP=" << std::scientific << std::setprecision(6)
		          << mc.aggregateImpactProb << " sigma=" << mc.aggregateImpactSigma << '\n';
		for (std::size_t e = 0; e < events.size(); ++e) {
			const double pe = mc.numConverged > 0 ? static_cast<double>(mc.hitsPerEvent[e]) / mc.numConverged : 0.0;
			std::cout << "  event[" << e << "]"
			          << " jd=" << std::fixed << std::setprecision(3) << events[e].jd << " hits=" << mc.hitsPerEvent[e]
			          << " p=" << std::scientific << std::setprecision(6) << pe << " minMissKm=" << std::fixed
			          << std::setprecision(1) << mc.minMissPerEventKm[e] << " medMissKm=" << mc.medianMissPerEventKm[e]
			          << '\n';
		}
	}

}  // namespace

[[nodiscard]] std::string describeMask(std::uint16_t mask)
{
	static constexpr const char* kNames[8] = {"Mercury", "Venus",  "Earth",  "Mars",
	                                          "Jupiter", "Saturn", "Uranus", "Neptune"};
	std::string out;
	const auto append = [&out](const char* name) {
		if (!out.empty())
			out += '+';
		out += name;
	};
	for (int bit = 0; bit < 8; ++bit) {
		if ((mask >> bit) & 1u)
			append(kNames[bit]);
	}
	if ((mask >> 11) & 1u)
		append("Moon");
	if ((mask >> 12) & 1u)
		append("Pluto");
	return out.empty() ? std::string("<none>") : out;
}

[[nodiscard]] const char* weightProfileName(orbital_core::AstrometricWeightProfile p)
{
	switch (p) {
		case orbital_core::AstrometricWeightProfile::Vfcc17:
			return "vfcc17";
		case orbital_core::AstrometricWeightProfile::Uniform:
			return "uniform";
	}
	return "unknown";
}

void printForceModel(const CliOptions& opt, const orbital_core::FitPipelineResult& result, bool sb441Available)
{
	const auto& req = opt.cfg.nbody;
	const auto& eff = result.effectiveConfig.nbody;
	const bool sb441Live = eff.useSb441 && sb441Available;

	std::cout << "[ForceModel] fit   planetMask=0x" << std::hex << eff.planetMask << std::dec << " ["
	          << describeMask(eff.planetMask) << "]";
	if (eff.planetMask != req.planetMask) {
		std::cout << "  <- auto-reduced from 0x" << std::hex << req.planetMask << std::dec;
	}
	std::cout << '\n';
	std::cout << "[ForceModel] prop  planetMask=0x" << std::hex << req.planetMask << std::dec << " ["
	          << describeMask(req.planetMask) << "]"
	          << "  (--ca-scan / --lov-impact / --obs-mc / --state-at)" << '\n';
	std::cout << "[ForceModel] nbody=" << (eff.enabled ? "on" : "off") << " step=" << eff.stepSize << "d accuracy=fixed"
	          << '\n';
	std::cout << "[ForceModel] sb441 requested=" << (req.useSb441 ? "on" : "off")
	          << " effective=" << (eff.useSb441 ? "on" : "off") << " kernel=" << (sb441Available ? "loaded" : "ABSENT")
	          << " -> asteroids=" << (sb441Live ? "SB441-N16 (SPK)" : "Kepler-3 fallback (Ceres/Pallas/Vesta)") << '\n';
	if (eff.targetNaifId)
		std::cout << "[ForceModel] target NAIF=" << eff.targetNaifId << " (self-perturbation excluded when present)\n";
	std::cout << "[Astrometry] weightProfile=" << weightProfileName(opt.cfg.weightProfile)
	          << " observerSigmas=" << (opt.cfg.obsWeighting.useObserverSigmas ? "on" : "off")
	          << " debias=" << (opt.debias ? "on" : "off") << " weights=" << (opt.cfg.useWeights ? "on" : "off")
	          << " sigmaClip=" << (opt.cfg.sigmaClip.enabled ? "on" : "off") << "@" << opt.cfg.sigmaClip.thresholdNSigma
	          << "sigma\n";
	std::cout << "[Astrometry] biasTable="
	          << (opt.debias ? result.observations.preparation->biasTable.string() : std::string("<none>")) << '\n';
	const NonGravModel& ng = result.fit.finalNonGrav;
	std::cout << "[ForceModel] nonGrav=" << (ng.enabled() ? "on" : "off") << '\n';

	if (result.fit.uncertainty.hasElementCovariance || result.fit.uncertainty.cometaryElementCovariance.has_value()) {
		const auto& sg = result.fit.uncertainty.derivedSigmas;
		constexpr double deg = Constants::RAD_TO_DEG;
		std::cout << "[Sigma] 1-sigma at the output epoch (a,q,Q in AU; angles in deg)" << '\n'
		          << std::scientific << std::setprecision(4) << "[Sigma]   a    " << sg.sigma_a << "   e    "
		          << sg.sigma_e << "   q    " << sg.sigma_q << "   Q    " << sg.sigma_Q << '\n'
		          << "[Sigma]   i    " << sg.sigma_i * deg << "   node " << sg.sigma_Om * deg << "   peri "
		          << sg.sigma_w * deg << "   M    " << sg.sigma_M * deg << '\n'
		          << "[Sigma]   n    " << sg.sigma_n * deg << "   P    " << sg.sigma_P << "   Tp   " << sg.sigma_Tp
		          << "   1/a  " << sg.sigma_inv_a << '\n'
		          << std::defaultfloat;
	}
}

void printStateAt(const CliOptions& opt, const orbital_core::FitPipelineResult& result)
{
	const Orbit_p& orbit = result.fit.fitEpochOrbit ? result.fit.fitEpochOrbit : result.fit.orbit;
	if (!orbit)
		throw std::runtime_error("--state-at requires a fitted orbit");

	NBodyPropagator prop = makeSurveyPropagator(opt.cfg, result.fit, IntegratorKind::RA15, opt.caCfg.encounterRadiusAu);
	prop.setTargetNaifId(result.effectiveConfig.nbody.targetNaifId);

	const Vec3D R0 = orbit->getR();
	const Vec3D V0 = orbit->getV() * Constants::GAUSS_K;
	const double epochJd = orbit->getEpoch();

	std::vector<double> targets;
	targets.reserve(opt.stateAtMjd.size());
	for (const double mjd : opt.stateAtMjd)
		targets.push_back(mjd + 2400000.5);

	const auto states = prop.propagateAllWithVelocity(R0, V0, epochJd, targets);

	const auto old = std::cout.flags();
	std::cout << std::uppercase << std::scientific << std::setprecision(16);
	std::cout << "[StateAt] frame=ECLM J2000 heliocentric  units=AU, AU/day"
	          << "  scale=TDB\n";
	std::cout << "[StateAt] epoch_mjd " << (epochJd - 2400000.5) << '\n';
	std::cout << "[StateAt] R0 " << R0.x << ' ' << R0.y << ' ' << R0.z << '\n';
	std::cout << "[StateAt] V0 " << V0.x << ' ' << V0.y << ' ' << V0.z << '\n';
	for (std::size_t k = 0; k < states.size() && k < targets.size(); ++k) {
		const auto& st = states[k];
		std::cout << "[StateAt] target_mjd " << opt.stateAtMjd[k] << '\n';
		std::cout << "[StateAt] R " << st.R.x << ' ' << st.R.y << ' ' << st.R.z << '\n';
		std::cout << "[StateAt] V " << st.V.x << ' ' << st.V.y << ' ' << st.V.z << '\n';
	}
	std::cout.flags(old);
}

int main(int argc, char* argv[])
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif

	try {
		CliOptions opt = parseArgs(argc, argv);

		if (opt.verbose) {
			Logger::setCallback([](LogLevel, const std::string& msg) { std::cout << msg << '\n'; });
		}

		getExecutableDir() = std::filesystem::path(argv[0]).parent_path();

		DE440Source de440;
		setDefaultEphemerisSource(&de440);
		SpkReader sb441;
		setDefaultAsteroidSource(&sb441);

		const auto result = runPipeline(opt);
		if (!opt.lovInputPath.empty())
			dumpLovInput(opt.lovInputPath, result.fit);
		printSummary(result);
		printForceModel(opt, result, sb441.isAvailable());

		std::filesystem::path solutionDirectory;
		{
			const auto& observations = result.observations;
			orbital_core::ExportContext exportContext;
			exportContext.objectName =
			    observations.objectName.empty() ? observations.source.stem().string() : observations.objectName;
			exportContext.fitResult = &result.fit;
			exportContext.observations = observations.all;
			exportContext.gaiaPhotocenterRadiusKm = observations.gaiaPhotocenterRadiusKm;
			exportContext.weightProfile = opt.cfg.weightProfile;
			exportContext.useWeights = opt.cfg.useWeights;
			exportContext.cometPhotometry = result.usesCometPhotometry();
			exportContext.photometry = result.effectiveConfig.photometry;
			const auto exportReport = orbital_core::Exporter::save(exportContext, {orbital_core::ExportProduct::All});
			solutionDirectory = exportReport.directory;
			std::cout << "Fit export written to: " << exportReport.directory.string() << '\n';
		}

		const bool needSurvey = opt.closeApproachScan || opt.lovImpact || opt.obsMcClones > 0;
		if (needSurvey) {
			Orbit_p surveyOrbit = result.fit.fitEpochOrbit ? result.fit.fitEpochOrbit : result.fit.orbit;
			if (!surveyOrbit)
				throw std::runtime_error("close-approach survey requires a fitted orbit");

			NBodyPropagator surveyProp =
			    makeSurveyPropagator(opt.cfg, result.fit, opt.surveyIntegrator, opt.caCfg.encounterRadiusAu);
			surveyProp.setTargetNaifId(result.effectiveConfig.nbody.targetNaifId);
			std::cout << "[SurveyIntegrator] name=" << surveyProp.integratorName()
			          << " accuracy=" << surveyProp.fixedAccuracyParameter() << '\n';

			std::vector<orbital_core::hazard::CloseApproachEvent> caEvents;
			if (opt.closeApproachScan || opt.obsMcClones > 0) {
				caEvents = orbital_core::hazard::findCloseApproaches(*surveyOrbit, surveyProp, opt.caCfg, nullptr);
			}

			if (opt.closeApproachScan)
				printCloseApproaches(caEvents);

			if (opt.lovImpact) {
				auto lovCfg = opt.lovCfg;
				lovCfg.closeApproach = opt.caCfg;
				printLovImpact(*surveyOrbit, result.fit, result.observations.all, result.effectiveConfig, surveyProp,
				               lovCfg);
			}

			if (opt.obsMcClones > 0) {
				printObsMc(*surveyOrbit, result.fit, surveyProp, caEvents, opt.obsMcClones, opt.obsMcSeed,
				           opt.obsMcAltitudeKm, opt.obsMcOutput, solutionDirectory, result.observations.objectName);
			}
		}

		if (!opt.stateAtMjd.empty())
			printStateAt(opt, result);

		if (!opt.residualsPath.empty()) {
			if (!result.fit.residuals.save(opt.residualsPath))
				throw std::runtime_error("failed to save residuals: " + opt.residualsPath.string());
			std::cout << "Residuals saved to: " << opt.residualsPath.string() << '\n';
		}

		return result.fit.orbit ? 0 : 2;
	} catch (const std::exception& e) {
		std::cerr << "TrajectaCLI: " << e.what() << '\n';
		return 1;
	}
}
