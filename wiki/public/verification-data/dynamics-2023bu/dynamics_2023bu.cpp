#include "NBodyPropagator.h"
#include "earth_orientation/EarthOrientation.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string_view>
#include <vector>

struct Reference {
	double mjd;
	Vec3D r, v;
};

int main(int argc, char** argv)
{
	if (argc < 5 || argc > 6 || (argc == 6 && std::string_view(argv[5]) != "--no-pluto")) {
		std::cerr << "Usage: dynamics_2023bu DE440 SB441 HORIZONS_2023 HORIZONS_100Y [--no-pluto]\n";
		return 2;
	}
	const bool usePluto = argc != 6;
	getExecutableDir() = ".";
	std::vector<Reference> refs;
	auto append = [&](const char* path) {
		std::ifstream input(path);
		const auto json = nlohmann::json::parse(input);
		const std::string text = json.at("result").get<std::string>();
		const auto start = text.find("$$SOE");
		const auto finish = text.find("$$EOE");
		if (start == std::string::npos || finish == std::string::npos)
			throw std::runtime_error("missing vectors");
		std::istringstream lines(text.substr(start + 5, finish - start - 5));
		std::string line, pos, vel;
		while (std::getline(lines, line)) {
			double jd;
			if (std::sscanf(line.c_str(), "%lf =", &jd) != 1 || jd < 2400000)
				continue;
			if (!std::getline(lines, pos) || !std::getline(lines, vel))
				throw std::runtime_error("missing vector line");
			Reference ref{};
			ref.mjd = jd - 2400000.5;
			if (std::sscanf(pos.c_str(), " X = %lf Y = %lf Z = %lf", &ref.r.x, &ref.r.y, &ref.r.z) != 3 ||
			    std::sscanf(vel.c_str(), " VX= %lf VY= %lf VZ= %lf", &ref.v.x, &ref.v.y, &ref.v.z) != 3)
				throw std::runtime_error("invalid vectors");
			refs.push_back(ref);
		}
	};
	for (int i = 3; i < 5; ++i)
		append(argv[i]);
	std::sort(refs.begin(), refs.end(), [](const auto& a, const auto& b) { return a.mjd < b.mjd; });
	refs.erase(std::unique(refs.begin(), refs.end(), [](const auto& a, const auto& b) { return a.mjd == b.mjd; }),
	           refs.end());
	if (refs.size() < 2) {
		std::cerr << "Unexpected reference count " << refs.size() << '\n';
		return 4;
	}
	DE440Source de(argv[1]);
	SpkReader sb(argv[2]);
	if (!de.isAvailable() || !sb.isAvailable())
		return 5;
	setDefaultEphemerisSource(&de);
	setDefaultAsteroidSource(&sb);
	const double tt = orbital_core::time::fromMjdTdb(refs[2].mjd, orbital_core::time::Scale::TT);
	const auto eop = orbital_core::earthOrientation().lookupMjdTt(tt);
	std::cout << "EOP available=" << orbital_core::earthOrientation().isAvailable()
	          << " records=" << orbital_core::earthOrientation().size() << " valid=" << eop.valid
	          << " extrapolated=" << eop.extrapolated << '\n';
	std::vector<double> times;
	for (const auto& ref : refs)
		times.push_back(ref.mjd);
	std::vector<NBodyPropagator::StateVel> first;
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		NBodyPropagator base(kind);
		base.setUseSb441(true);
		base.setPlanetEnabled(EphemerisPlanet::Pluto, usePluto);
		auto p = base;
		const auto out = p.propagateAllWithVelocityMjd(refs[0].r, refs[0].v, refs[0].mjd, times, p.nonGrav());
		if (out.size() != refs.size())
			return 6;
		std::cout << "INTEGRATOR " << p.integratorName() << " NG=" << p.nonGrav().enabled() << " PLUTO=" << usePluto
		          << '\n';
		for (std::size_t i = 0; i < out.size(); ++i) {
			const Vec3D dr = (out[i].R - refs[i].r) * 149597870700.0;
			const Vec3D dv = (out[i].V - refs[i].v) * (149597870700.0 / 86400.0);
			std::cout << std::setprecision(14) << "RESULT " << p.integratorName() << ' ' << refs[i].mjd << " dr_m "
			          << dr.magnitude() << " dv_m_s " << dv.magnitude() << " xyz_m " << dr.x << ' ' << dr.y << ' '
			          << dr.z << " vxyz_m_s " << dv.x << ' ' << dv.y << ' ' << dv.z;
			if (!first.empty())
				std::cout << " integrator_diff_m " << (out[i].R - first[i].R).magnitude() * 149597870700.0;
			std::cout << '\n';
		}
		if (first.empty())
			first = out;
	}
}
