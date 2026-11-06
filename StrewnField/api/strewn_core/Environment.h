#pragma once

#include "strewn_core/Export.h"
#include "strewn_core/Types.h"

namespace strewn {

	struct GeodeticPosition {
		double latitudeRad{};
		double longitudeRad{};
		double heightWgs84M{};
	};

	struct AirSample {
		double densityKgM3{};
		double temperatureK{};
		Vector3 windEnuMps{};
	};

	struct AtmosphereGrid {
		std::vector<double> latitudesDeg{0.0};
		std::vector<double> longitudesDeg{0.0};
		double firstHeightM{};
		double heightStepM{50.0};
		std::size_t heightCount{};
		// Latitude, longitude, height; singleton horizontal axes describe a column.
		std::vector<AirSample> samples;
	};

	struct TerrainGrid {
		std::vector<double> latitudesDeg{0.0};
		std::vector<double> longitudesDeg{0.0};
		std::vector<double> heightsWgs84M{0.0};
	};

	STREWN_CORE_EXPORT GeodeticPosition geodetic(const Vector3& ecefM);
	STREWN_CORE_EXPORT Vector3 ecefPosition(const GeodeticPosition& position);
	STREWN_CORE_EXPORT Vector3 enuToEcef(const Vector3& enu, double latitudeRad, double longitudeRad);
	STREWN_CORE_EXPORT AirSample standardAtmosphere(double heightM);

	class STREWN_CORE_EXPORT Environment {
	   public:
		Environment(AtmosphereGrid atmosphere, TerrainGrid terrain);
		static Environment idealized(double groundHeightWgs84M = 0.0);
		AirSample air(const GeodeticPosition& position) const;
		double ground(const GeodeticPosition& position) const;
		double clearance(const GeodeticPosition& position) const;

		const AtmosphereGrid& atmosphereGrid() const noexcept { return m_atmosphere; }

		const TerrainGrid& terrainGrid() const noexcept { return m_terrain; }

	   private:
		AtmosphereGrid m_atmosphere;
		TerrainGrid m_terrain;
		double m_maxGroundHeightM{};
	};

}  // namespace strewn
