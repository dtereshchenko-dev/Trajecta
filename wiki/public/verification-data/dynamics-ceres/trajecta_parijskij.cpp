#include "NBodyPropagator.h"
#include "integrators/IAS15Integrator.h"
#include "integrators/RA15Integrator.h"
#include "integrators/BarycentricIntegrator.h"
#include "propagator/NBodyDynamics.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
	if (argc < 9 || argc > 12) {
		std::cerr
		    << "usage: trajecta_parijskij DE440 SB441 x y z vx vy vz "
		       "[--exclude-ceres | --adaptive-final | --ias15-final | --ias15-grid | --ias15-trace | --ias15-steps | "
		       "--fixed-step days [--sample-grid] | --ias15-fixed-step days [--sample-grid] | "
		       "--ias15-initial-step days | --ra15-initial-step days | --ra15-trace]\n";
		return 2;
	}
	DE440Source de(argv[1]);
	SpkReader sb(argv[2]);
	if (!de.isAvailable() || !sb.isAvailable()) {
		std::cerr << "ephemeris unavailable\n";
		return 3;
	}
	setDefaultEphemerisSource(&de);
	setDefaultAsteroidSource(&sb);
	const Vec3D r(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]));
	const Vec3D v(std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]));
	constexpr double epoch = 49718.0;
	constexpr int samples = 3653;
	const std::string mode = argc >= 10 ? argv[9] : "";
	const bool customInitial = argc == 11 && mode == "--ias15-initial-step";
	const bool raAdaptive = (argc == 11 && mode == "--ra15-initial-step") || (argc == 10 && mode == "--ra15-trace");
	const bool iasAdaptive = customInitial || (argc == 10 && (mode == "--ias15-final" || mode == "--ias15-grid" ||
	                                                          mode == "--ias15-trace" || mode == "--ias15-steps"));
	if (argc >= 11 || iasAdaptive || raAdaptive) {
		const bool fixed = !iasAdaptive && !raAdaptive;
		const bool ias = iasAdaptive || mode == "--ias15-fixed-step";
		if (fixed && mode != "--fixed-step" && mode != "--ias15-fixed-step")
			return 2;
		if (argc == 12 && std::string(argv[11]) != "--sample-grid")
			return 2;
		const double step = (fixed || customInitial || (raAdaptive && argc == 11)) ? std::stod(argv[10]) : 10.0;
		if (!(step > 0.0) || !std::isfinite(step))
			return 2;
		NBodyDynamics dynamics;
		dynamics.setUseSb441(true);
		std::unique_ptr<IIntegrator> engine = ias ? std::unique_ptr<IIntegrator>(std::make_unique<IAS15Integrator>())
		                                          : std::unique_ptr<IIntegrator>(std::make_unique<RA15Integrator>());
		BarycentricIntegrator integrator(std::move(engine));
		IntegrationOptions options;
		options.initialStepDays = step;
		if (fixed)
			options.nextEncounterStep = [step](double, double, const IntegratedState&) {
				return std::optional<double>{step};
			};
		if (mode == "--ias15-steps")
			options.nextEncounterStepOffset = [](double step, double, double offset, const IntegratedState&) {
				std::cerr << std::setprecision(17) << "step " << offset << ' ' << step << '\n';
				return std::optional<double>{};
			};
		std::vector<double> target;
		if (argc == 12 || mode == "--ias15-grid" || mode == "--ias15-trace" || mode == "--ra15-trace") {
			target.reserve(samples);
			for (int i = 0; i < samples; ++i)
				target.push_back(epoch + 3653.0 * i / (samples - 1));
		} else {
			target.push_back(epoch + 3653.0);
		}
		std::vector<IntegratedState> out(target.size());
		const auto report = integrator.integrateState(dynamics, NonGravModel{}, {r, v}, epoch, target, out, options);
		std::cerr << "accepted_steps=" << report.acceptedSteps << " rejected_steps=" << report.rejectedSteps
		          << " unconverged_reported=" << report.unconvergedEncounterSteps << '\n';
		if (report.completedOutputs != target.size() || report.unconvergedEncounterSteps != 0)
			return 5;
		std::cout << std::setprecision(17);
		const std::size_t first = (mode == "--ias15-trace" || mode == "--ra15-trace") ? 0 : out.size() - 1;
		for (std::size_t i = first; i < out.size(); ++i) {
			const auto& state = out[i];
			std::cout << target[i] - epoch << ' ' << state.R.x << ' ' << state.R.y << ' ' << state.R.z << ' '
			          << state.V.x << ' ' << state.V.y << ' ' << state.V.z << '\n';
		}
		return 0;
	}
	std::vector<double> targets;
	targets.reserve(samples);
	for (int i = 0; i < samples; ++i)
		targets.push_back(epoch + 3653.0 * i / (samples - 1));

	NBodyPropagator prop(IntegratorKind::RA15);
	prop.setUseSb441(true);
	prop.setPlanetEnabled(EphemerisPlanet::Pluto, true);
	if (argc == 10) {
		if (mode == "--exclude-ceres")
			prop.setTargetNaifId(2000001);
		else if (mode == "--adaptive-final")
			targets = {epoch + 3653.0};
		else
			return 2;
	}
	const auto rows = prop.propagateAllWithVelocityMjd(r, v, epoch, targets, prop.nonGrav());
	if (rows.size() != targets.size()) {
		std::cerr << "incomplete propagation: " << rows.size() << " of " << targets.size() << '\n';
		return 4;
	}
	std::cout << std::setprecision(17);
	for (std::size_t i = 0; i < rows.size(); ++i) {
		const auto& s = rows[i];
		std::cout << targets[i] - epoch << ' ' << s.R.x << ' ' << s.R.y << ' ' << s.R.z << ' ' << s.V.x << ' ' << s.V.y
		          << ' ' << s.V.z << '\n';
	}
}
