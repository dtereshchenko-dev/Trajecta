#pragma once

#include <QQmlEngine>
#include <QVector3D>
#include <orbital_core/Orbit.h>
#include "util/PlanetPositionsProvider.h"

class OrbitSceneModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the solution")
	Q_PROPERTY(QVariantList planets READ planets NOTIFY changed)
	Q_PROPERTY(QVariantList planetPositions READ planetPositions NOTIFY positionChanged)
	Q_PROPERTY(QVector3D objectPosition READ objectPosition NOTIFY positionChanged)
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(bool canAnimate READ canAnimate NOTIFY changed)
	Q_PROPERTY(double frameRadius READ frameRadius NOTIFY changed)
	Q_PROPERTY(double epochMjdTdb READ epochMjdTdb NOTIFY positionChanged)
	Q_PROPERTY(QString epochUtc READ epochUtc NOTIFY positionChanged)
	Q_PROPERTY(double position READ position WRITE setPosition NOTIFY positionChanged)
	Q_PROPERTY(double spanDays READ spanDays NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY positionChanged)

   public:
	explicit OrbitSceneModel(QObject* parent = nullptr);
	void setOrbit(const Orbit_p& orbit);

	QVariantList planets() const { return m_planets; }

	QVariantList planetPositions() const { return m_planetPositions; }

	QVector3D objectPosition() const { return m_objectPosition; }

	bool available() const { return bool(m_orbit); }

	bool canAnimate() const { return m_canAnimate; }

	double frameRadius() const { return m_frameRadius; }

	double epochMjdTdb() const;
	QString epochUtc() const;

	double position() const { return m_position; }

	double spanDays() const { return m_spanDays; }

	QString error() const { return m_error; }

	void setPosition(double position);

   signals:
	void changed();
	void positionChanged();

   private:
	void updatePlanets();
	Orbit_p m_orbit;
	PlanetPositionsProvider m_planetProvider;
	QVariantList m_planets;
	QVariantList m_planetPositions;
	QVector3D m_objectPosition;
	double m_position{};
	double m_spanDays{};
	double m_frameRadius{2};
	bool m_canAnimate{};
	QString m_error;
};
