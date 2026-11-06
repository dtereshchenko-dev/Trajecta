#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/Constants.h"

namespace orbital_core::hazard {
	namespace {


		constexpr PlanetData kMercury{22'031.868551400003, 2'440.53, 0.0, 60.0e-6};
		constexpr PlanetData kVenus{324'858.592, 6'051.8, 0.0, 4.458e-6};
		constexpr PlanetData kEarth{398'600.43550702254, Constants::EARTH_EQ_RADIUS_KM, Constants::EARTH_FLATTENING,
		                            1.08263e-3};
		constexpr PlanetData kMars{42'828.375815756102, 3'396.19, 0.0, 1.96045e-3};
		constexpr PlanetData kJupiter{126'712'764.09999998, 71'492.0, 0.0, 1.4736e-2};
		constexpr PlanetData kSaturn{37'940'584.841799997, 60'268.0, 0.0, 1.6298e-2};
		constexpr PlanetData kUranus{5'794'556.3999999985, 25'559.0, 0.0, 3.343e-3};
		constexpr PlanetData kNeptune{6'836'527.100580399, 24'764.0, 0.0, 3.411e-3};
		constexpr PlanetData kMoon{4'902.80011845755, 1'737.4, 0.0, 2.027e-4};
	}  // namespace

	const PlanetData& getPlanetData(Planet p) noexcept
	{
		switch (p) {
			case Planet::Mercury:
				return kMercury;
			case Planet::Venus:
				return kVenus;
			case Planet::Earth:
				return kEarth;
			case Planet::Mars:
				return kMars;
			case Planet::Jupiter:
				return kJupiter;
			case Planet::Saturn:
				return kSaturn;
			case Planet::Uranus:
				return kUranus;
			case Planet::Neptune:
				return kNeptune;
			case Planet::Moon:
				return kMoon;
			case Planet::Count_:
				break;
		}
		return kEarth;
	}
}  // namespace orbital_core::hazard
