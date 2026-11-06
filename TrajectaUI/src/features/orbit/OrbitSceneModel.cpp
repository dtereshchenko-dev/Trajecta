#include "OrbitSceneModel.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include "util/FormatUtils.h"
#include "util/CalculationTime.h"

OrbitSceneModel::OrbitSceneModel(QObject* parent) : QObject(parent) {}

void OrbitSceneModel::setOrbit(const Orbit_p& orbit)
{
	m_orbit = orbit ? std::make_shared<Orbit>(*orbit) : nullptr;
	m_position = 0;
	m_spanDays = 0;
	m_canAnimate = false;
	m_frameRadius = 2;
	m_objectPosition = {};
	m_planets.clear();
	m_planetPositions.clear();
	m_error.clear();
	if (m_orbit) {
		const double e = m_orbit->getExcentricity();
		m_objectPosition = FormatUtils::eclToQt3D(m_orbit->getR());
		m_canAnimate = std::isfinite(m_orbit->getSemiMajorAxis()) && std::isfinite(m_orbit->getMeanAnomaly()) &&
		               std::isfinite(e) && std::abs(e - 1) > 1e-5 && m_orbit->getPerigee() > 0;
		const double period = m_orbit->getOrbitalPeriod();
		m_spanDays = m_canAnimate ? (std::isfinite(period) && period > 0 ? std::min(period, 36525.0) : 365.25) : 0;
		const double extent = e < 0.99 ? m_orbit->getApogee() : m_orbit->getR().magnitude() * 1.5;
		if (std::isfinite(extent))
			m_frameRadius = std::clamp(extent, 1.6, 100000.0);
		updatePlanets();
		m_planets = m_planetProvider.positions();
	}
	emit changed();
	emit positionChanged();
}

double OrbitSceneModel::epochMjdTdb() const { return m_orbit ? m_orbit->epochMjd() + m_position * m_spanDays : 0; }

QString OrbitSceneModel::epochUtc() const
{
	if (!m_orbit)
		return {};
	return CalculationTime::formatUtc(orbital_core::time::fromMjdTdb(epochMjdTdb(), orbital_core::time::Scale::UTC));
}

void OrbitSceneModel::setPosition(double position)
{
	if (!m_orbit || !m_canAnimate || !std::isfinite(position))
		return;
	position = std::clamp(position, 0.0, 1.0);
	if (position == m_position)
		return;
	const auto r = position == 0
	                   ? m_orbit->getR()
	                   : Calculations::positionAt(m_orbit->getSemiMajorAxis(), m_orbit->getExcentricity(),
	                                              m_orbit->getInclination(), m_orbit->getLongitudeOfAscendingNode(),
	                                              m_orbit->getArgumentOfPericentre(), m_orbit->getMeanAnomaly(),
	                                              m_orbit->epochMjd(), m_orbit->epochMjd() + position * m_spanDays,
	                                              Constants::GAUSS_K * Constants::GAUSS_K, m_orbit->getType());
	if (!std::isfinite(r.x) || !std::isfinite(r.y) || !std::isfinite(r.z)) {
		m_error = tr("Two-body position unavailable");
		emit positionChanged();
		return;
	}
	m_position = position;
	m_objectPosition = FormatUtils::eclToQt3D(r);
	updatePlanets();
	emit positionChanged();
}

void OrbitSceneModel::updatePlanets()
{
	m_planetPositions.clear();
	m_error.clear();
	try {
		m_planetProvider.setEpoch(epochMjdTdb() + 2400000.5);
		const auto positions = m_planetProvider.positions();
		for (const auto& entry : positions)
			m_planetPositions.append(entry.toMap().value("pos"));
		if (m_planetPositions.isEmpty())
			m_error = tr("Planet ephemerides unavailable");
	} catch (const std::exception&) {
		m_error = tr("Planet ephemerides unavailable at this epoch");
	}
}
