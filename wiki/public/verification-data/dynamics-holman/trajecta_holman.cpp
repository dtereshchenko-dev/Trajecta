#include "NBodyPropagator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
	if (argc != 9 && argc != 10) {
		std::cerr << "usage: trajecta_holman DE440 SB441 x y z vx vy vz [--ias15]\n";
		return 2;
	}
	const bool ias = argc == 10 && std::string(argv[9]) == "--ias15";
	if (argc == 10 && !ias)
		return 2;
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
	constexpr double epoch = 58849.0;
	const std::vector<double> days{0, 1, 100, 1000, 5000, 10000};
	std::vector<double> targets;
	for (double day : days)
		targets.push_back(epoch + day);
	NBodyPropagator prop(ias ? IntegratorKind::IAS15 : IntegratorKind::RA15);
	prop.setUseSb441(true);
	prop.setPlanetEnabled(EphemerisPlanet::Pluto, true);
	IntegrationReport report;
	const auto rows = prop.propagateAllWithStmMjd(r, v, epoch, targets, prop.nonGrav(), &report);
	if (rows.size() != days.size()) {
		std::cerr << "incomplete propagation: " << rows.size() << " of " << days.size() << ", reached MJD "
		          << report.reachedMjd << '\n';
		return 4;
	}
	std::cout << std::setprecision(17);
	for (std::size_t i = 0; i < rows.size(); ++i) {
		const auto& s = rows[i];
		std::cout << days[i] << ' ' << s.R.x << ' ' << s.R.y << ' ' << s.R.z << ' ' << s.V.x << ' ' << s.V.y << ' '
		          << s.V.z;
		for (int row = 0; row < 6; ++row)
			for (int col = 0; col < 6; ++col)
				std::cout << ' ' << s.Phi(row, col);
		std::cout << '\n';
	}
}
