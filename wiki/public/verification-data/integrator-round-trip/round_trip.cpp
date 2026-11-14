#include "NBodyPropagator.h"
#include "earth_orientation/EarthOrientation.h"
#include "integrators/IAS15Integrator.h"
#include "integrators/RA15Integrator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "orbital_core/forces/SunGravity.h"

#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
	constexpr double kAuMeters = 149597870700.0;
	constexpr double kSecondsPerDay = 86400.0;

	class KeplerDynamics final : public IIntegrationDynamics {
	   public:
		Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D&, double, const NonGravModel&) const override
		{
			return r * (-Constants::GAUSS_K * Constants::GAUSS_K / std::pow(r.magnitude(), 3));
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D&, double, const NonGravModel&, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			orbital_core::forces::ForceContext c{};
			c.R = r;
			c.centralNewtonianGm = Constants::GAUSS_K * Constants::GAUSS_K;
			a = {};
			dr.setZero();
			dv.setZero();
			orbital_core::forces::SunGravity{}.evaluate(c, a, dr, dv);
		}
	};

	struct Initial {
		double mjd{};
		Vec3D r, v;
	};

	Initial readInitial(const char* path)
	{
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error("Cannot open reference file");
		const auto json = nlohmann::json::parse(input);
		const std::string text = json.at("result").get<std::string>();
		const auto start = text.find("$$SOE");
		const auto end = text.find("$$EOE");
		if (start == std::string::npos || end == std::string::npos || end <= start)
			throw std::runtime_error("Missing reference vectors");
		std::istringstream lines(text.substr(start + 5, end - start - 5));
		std::string line, pos, vel;
		while (std::getline(lines, line)) {
			double jd;
			if (std::sscanf(line.c_str(), "%lf =", &jd) != 1 || jd < 2400000)
				continue;
			if (!std::getline(lines, pos) || !std::getline(lines, vel))
				throw std::runtime_error("Incomplete vector");
			Initial result;
			result.mjd = jd - 2400000.5;
			if (std::sscanf(pos.c_str(), " X = %lf Y = %lf Z = %lf", &result.r.x, &result.r.y, &result.r.z) != 3 ||
			    std::sscanf(vel.c_str(), " VX= %lf VY= %lf VZ= %lf", &result.v.x, &result.v.y, &result.v.z) != 3)
				throw std::runtime_error("Invalid vector");
			return result;
		}
		throw std::runtime_error("No initial vector");
	}

	void printResult(const char* scenario, const char* method, double days, const Vec3D& initialR,
	                 const Vec3D& initialV, const Vec3D& returnedR, const Vec3D& returnedV)
	{
		std::cout << std::setprecision(17) << "RESULT " << scenario << ' ' << method << ' ' << days << " position_m "
		          << (returnedR - initialR).magnitude() * kAuMeters << " velocity_m_s "
		          << (returnedV - initialV).magnitude() * kAuMeters / kSecondsPerDay << '\n';
	}

	void runKepler(IntegratorKind kind, const char* name, double speedFactor, double days, bool circularReference)
	{
		auto makeEngine = [&]() -> std::unique_ptr<IIntegrator> {
			if (kind == IntegratorKind::IAS15)
				return std::make_unique<IAS15Integrator>();
			return std::make_unique<RA15Integrator>();
		};
		const double epoch = 60000.0;
		const IntegratedState initial{Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, Constants::GAUSS_K * speedFactor, 0.0)};
		KeplerDynamics dynamics;
		NonGravModel ng;
		ng.setEnabled(false);
		IntegrationOptions options;
		options.initialStepDays = 10.0;
		auto propagate = [&](const IntegratedState& state, double from, double to) {
			auto engine = makeEngine();
			std::array<double, 1> targets{to};
			std::array<IntegratedState, 1> out{};
			const auto report = engine->integrateState(dynamics, ng, state, from, targets, out, options);
			if (report.completedOutputs != 1)
				throw std::runtime_error("Kepler propagation failed");
			return out[0];
		};
		const auto forward = propagate(initial, epoch, epoch + days);
		const auto returned = propagate(forward, epoch + days, epoch);
		printResult(circularReference ? "circular" : "eccentric", name, days, initial.R, initial.V, returned.R,
		            returned.V);
		if (circularReference) {
			const double angle = Constants::GAUSS_K * ((epoch + days) - epoch);
			const Vec3D r(std::cos(angle), std::sin(angle), 0.0);
			const Vec3D v(-Constants::GAUSS_K * std::sin(angle), Constants::GAUSS_K * std::cos(angle), 0.0);
			printResult("circular_analytic", name, days, r, v, forward.R, forward.V);
		}
	}

	void runNBody(IntegratorKind kind, const char* name, const char* object, const Initial& initial, double days,
	              bool apophis)
	{
		auto setup = [&](NBodyPropagator& p) {
			p.setUseSb441(true);
			p.setPlanetEnabled(EphemerisPlanet::Pluto, true);
			if (apophis) {
				p.nonGrav().setModel(NgModel::InverseSquare);
				p.nonGrav().setEnabled(true);
				p.nonGrav().setParams(4.999999873689e-13, -2.901766720242e-14, 0.0);
			} else {
				p.nonGrav().setEnabled(false);
			}
		};
		NBodyPropagator forward(kind), backward(kind);
		setup(forward);
		setup(backward);
		const auto propagated = forward.propagateAllWithVelocityMjd(initial.r, initial.v, initial.mjd,
		                                                            {initial.mjd + days}, forward.nonGrav());
		if (propagated.size() != 1)
			throw std::runtime_error("Forward N-body propagation failed");
		const auto returned = backward.propagateAllWithVelocityMjd(propagated[0].R, propagated[0].V, initial.mjd + days,
		                                                           {initial.mjd}, backward.nonGrav());
		if (returned.size() != 1)
			throw std::runtime_error("Backward N-body propagation failed");
		printResult(object, name, days, initial.r, initial.v, returned[0].R, returned[0].V);
	}
}  // namespace

int main(int argc, char** argv)
{
	if (argc != 5) {
		std::cerr << "Usage: round_trip DE440 SB441 HORIZONS_BU HORIZONS_APOPHIS\n";
		return 2;
	}
	try {
		const double period = 2.0 * Constants::PI / Constants::GAUSS_K;
		for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
			const char* name = kind == IntegratorKind::RA15 ? "RA15" : "IAS15";
			runKepler(kind, name, 1.0, period, true);
			runKepler(kind, name, 1.303840, 100.0 * period, false);
		}
		getExecutableDir() = ".";
		DE440Source de(argv[1]);
		SpkReader sb(argv[2]);
		if (!de.isAvailable() || !sb.isAvailable())
			throw std::runtime_error("Ephemeris unavailable");
		setDefaultEphemerisSource(&de);
		setDefaultAsteroidSource(&sb);
		const auto bu = readInitial(argv[3]);
		const auto apophis = readInitial(argv[4]);
		std::cout << "EOP available=" << orbital_core::earthOrientation().isAvailable()
		          << " records=" << orbital_core::earthOrientation().size() << '\n';
		for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
			const char* name = kind == IntegratorKind::RA15 ? "RA15" : "IAS15";
			for (double days : {7.0, 8.0, 40.0, 365.0, 1000.0})
				runNBody(kind, name, "2023BU", bu, days, false);
			for (double days : {102.0, 103.0, 110.0, 365.0, 1000.0})
				runNBody(kind, name, "Apophis", apophis, days, true);
		}
		setDefaultAsteroidSource(nullptr);
		setDefaultEphemerisSource(nullptr);
	} catch (const std::exception& e) {
		std::cerr << e.what() << '\n';
		return 1;
	}
}
