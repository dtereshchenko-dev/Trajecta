#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QSettings>
#include <orbital_core/FitConfig.h>

class PhotometrySettings : public QObject {
	Q_OBJECT
	QML_ELEMENT
	Q_PROPERTY(double asteroidG MEMBER m_asteroidG NOTIFY changed)
	Q_PROPERTY(double cometK1 MEMBER m_cometK1 NOTIFY changed)
	Q_PROPERTY(double cometK2 MEMBER m_cometK2 NOTIFY changed)
	Q_PROPERTY(double nuclearPhaseCoefficient MEMBER m_nuclearPhaseCoefficient NOTIFY changed)

   public:
	using QObject::QObject;

	[[nodiscard]] orbital_core::PhotometryOptions options() const noexcept
	{
		return {m_asteroidG, m_cometK1, m_cometK2, m_nuclearPhaseCoefficient};
	}

	void save(QSettings& settings) const;
	void load(QSettings& settings);

   signals:
	void changed();

   private:
	double m_asteroidG{orbital_core::PhotometryOptions{}.asteroidG};
	double m_cometK1{orbital_core::PhotometryOptions{}.cometK1};
	double m_cometK2{orbital_core::PhotometryOptions{}.cometK2};
	double m_nuclearPhaseCoefficient{orbital_core::PhotometryOptions{}.cometNuclearPhaseCoefficient};
};
