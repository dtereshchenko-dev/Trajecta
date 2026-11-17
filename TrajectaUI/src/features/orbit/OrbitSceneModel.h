#pragma once

#include <QQmlEngine>
#include <QVector3D>
#include <orbital_core/Orbit.h>
#include "RelativeMotionModel.h"
#include "util/PlanetPositionsProvider.h"

namespace orbital_core {
	struct FitPipelineResult;
}

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
	Q_PROPERTY(QString epochDisplay READ epochDisplay NOTIFY positionChanged)
	Q_PROPERTY(QString epochDisplayScale READ epochDisplayScale NOTIFY changed)
	Q_PROPERTY(double position READ position WRITE setPosition NOTIFY positionChanged)
	Q_PROPERTY(double spanDays READ spanDays NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY positionChanged)
	Q_PROPERTY(bool relativeMode READ relativeMode WRITE setRelativeMode NOTIFY changed)
	Q_PROPERTY(RelativeMotionModel* relative READ relative CONSTANT)

   public:
	explicit OrbitSceneModel(QObject* parent = nullptr);
	void setOrbit(const Orbit_p& orbit);
	void setRelativeSolution(const orbital_core::FitPipelineResult& result);

	RelativeMotionModel* relative() { return &m_relative; }

	QVariantList planets() const { return m_planets; }

	QVariantList planetPositions() const { return m_planetPositions; }

	QVector3D objectPosition() const;

	bool available() const { return bool(m_orbit); }

	bool canAnimate() const { return m_relativeMode ? m_relative.relativeAvailable() : m_canAnimate; }

	double frameRadius() const;
	double epochMjdTdb() const;
	QString epochUtc() const;
	QString epochDisplay() const;
	QString epochDisplayScale() const;

	double position() const { return m_position; }

	double spanDays() const;

	QString error() const { return m_error; }

	void setPosition(double position);

	bool relativeMode() const { return m_relativeMode; }

	void setRelativeMode(bool value);
	Q_INVOKABLE bool seekCaptureStart(int index);

   signals:
	void changed();
	void positionChanged();

   private:
	void updatePlanets();
	Orbit_p m_orbit;
	PlanetPositionsProvider m_planetProvider;
	RelativeMotionModel m_relative;
	QVariantList m_planets, m_planetPositions;
	QVector3D m_objectPosition;
	double m_position{}, m_spanDays{}, m_frameRadius{2};
	bool m_canAnimate{}, m_relativeMode{};
	QString m_error;
};
