#include "util/PlanetPositionsProvider.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Orbit.h"
#include "util/FormatUtils.h"
#include "orbital_core/Types.h"

void PlanetPositionsProvider::setEpoch(double jd)
{
	m_positions.clear();

	auto* src = dynamic_cast<DE440Source*>(getDefaultEphemerisSource());
	if (!src || !src->isAvailable()) {
		emit positionsChanged();
		return;
	}

	struct P {
		EphemerisPlanet id;
		const char* name;
		const char* color;
	};

	static constexpr P kPlanets[] = {
	    {EphemerisPlanet::Mercury, "Mercury", "#a0a0a0"}, {EphemerisPlanet::Venus, "Venus", "#e8c46a"},
	    {EphemerisPlanet::Earth, "Earth", "#4a7dff"},     {EphemerisPlanet::Mars, "Mars", "#e87461"},
	    {EphemerisPlanet::Jupiter, "Jupiter", "#d4a373"}, {EphemerisPlanet::Saturn, "Saturn", "#f0d080"},
	    {EphemerisPlanet::Uranus, "Uranus", "#7ec8e3"},   {EphemerisPlanet::Neptune, "Neptune", "#5b7dd4"},
	};

	for (const auto& [id, name, color] : kPlanets) {
		auto ephData = src->queryJD(id, jd);
		Vec3D pos = ephData.R;
		if (pos.magnitude() < 0.01)
			continue;

		QVariantMap entry{
		    {"name", name},
		    {"color", color},
		    {"pos", FormatUtils::eclToQt3D(pos)},
		};

		Vec3D vGauss = ephData.V / Constants::GAUSS_K;
		auto orbit = make_orbit_from_states_vectors(StateVectors(pos, vGauss), jd);
		if (orbit) {
			entry["a"] = orbit->getSemiMajorAxis();
			entry["e"] = orbit->getExcentricity();
			entry["i"] = Convertors::radiansToDegrees(orbit->getInclination());
			entry["Om"] = Convertors::radiansToDegrees(orbit->getLongitudeOfAscendingNode());
			entry["w"] = Convertors::radiansToDegrees(orbit->getArgumentOfPericentre());
		}

		m_positions.append(entry);
	}

	emit positionsChanged();
}
