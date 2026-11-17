#include "OrbitSceneModel.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <orbital_core/FitPipeline.h>
#include "util/FormatUtils.h"
#include "util/CalculationTime.h"

OrbitSceneModel::OrbitSceneModel(QObject* parent) : QObject(parent), m_relative(this)
{
	connect(&m_relative, &RelativeMotionModel::positionChanged, this, &OrbitSceneModel::positionChanged);
	connect(&m_relative, &RelativeMotionModel::resultReady, this, [this](bool accepted) {
		if (accepted)
			m_position = 0;
		emit changed();
	});
}

void OrbitSceneModel::setOrbit(const Orbit_p& orbit)
{
	m_relative.reset();
	m_relativeMode = false;
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

void OrbitSceneModel::setRelativeSolution(const orbital_core::FitPipelineResult& result)
{
	m_relative.setSolution(result);
}

QVector3D OrbitSceneModel::objectPosition() const
{
	return m_relativeMode && m_relative.relativeAvailable() ? m_relative.objectPosition() : m_objectPosition;
}

double OrbitSceneModel::frameRadius() const
{
	return m_relativeMode && m_relative.relativeAvailable() ? m_relative.frameRadius() : m_frameRadius;
}

double OrbitSceneModel::spanDays() const
{
	return m_relativeMode && m_relative.relativeAvailable() ? m_relative.spanDays() : m_spanDays;
}

double OrbitSceneModel::epochMjdTdb() const
{
	if (m_relativeMode && m_relative.relativeAvailable())
		return m_relative.firstMjdTdb() + m_position * spanDays();
	return m_orbit ? m_orbit->epochMjd() + m_position * m_spanDays : 0;
}

QString OrbitSceneModel::epochUtc() const
{
	if (!m_orbit)
		return {};
	return CalculationTime::formatUtc(orbital_core::time::fromMjdTdb(epochMjdTdb(), orbital_core::time::Scale::UTC));
}

QString OrbitSceneModel::epochDisplayScale() const
{
	return m_relativeMode && m_relative.relativeAvailable() ? m_relative.relativeResultTimeScale()
	                                                        : QStringLiteral("UTC");
}

QString OrbitSceneModel::epochDisplay() const
{
	if (!m_orbit)
		return {};
	const auto scale =
	    epochDisplayScale() == QLatin1String("TDB") ? orbital_core::time::Scale::TDB : orbital_core::time::Scale::UTC;
	return CalculationTime::formatUtc(orbital_core::time::fromMjdTdb(epochMjdTdb(), scale));
}

void OrbitSceneModel::setPosition(double position)
{
	if (!m_orbit || !canAnimate() || !std::isfinite(position))
		return;
	position = std::clamp(position, 0.0, 1.0);
	if (position == m_position)
		return;
	if (m_relativeMode) {
		m_position = position;
		m_relative.setPosition(position);
		return;
	}
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

bool OrbitSceneModel::seekCaptureStart(int index)
{
	if (!m_relativeMode || !m_relative.relativeAvailable())
		return false;
	const auto epoch = m_relative.captureStartEpoch(index);
	const double start = m_relative.firstMjdTdb();
	const double duration = spanDays();
	if (!epoch || !std::isfinite(*epoch) || !std::isfinite(duration) || duration <= 0)
		return false;
	setPosition((*epoch - start) / duration);
	return true;
}

void OrbitSceneModel::setRelativeMode(bool value)
{
	if (m_relativeMode == value)
		return;
	m_relativeMode = value;
	m_position = 0;
	if (value && m_relative.relativeAvailable())
		m_relative.setPosition(0);
	else if (!value && m_orbit) {
		m_objectPosition = FormatUtils::eclToQt3D(m_orbit->getR());
		updatePlanets();
	}
	emit changed();
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
