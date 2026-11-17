#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QSettings>

#include "NBodyPropagator.h"

class PerturbationsSettings : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(bool enabled MEMBER m_enabled NOTIFY changed)
	Q_PROPERTY(uint planetMask MEMBER m_planetMask NOTIFY changed)
	Q_PROPERTY(bool useSb441 MEMBER m_useSb441 NOTIFY changed)
	Q_PROPERTY(double integratorInitialStepDays MEMBER m_integratorInitialStepDays NOTIFY changed)
	Q_PROPERTY(int relativityMode MEMBER m_relativityMode NOTIFY changed)
	Q_PROPERTY(bool earthJ2 MEMBER m_earthJ2 NOTIFY changed)
	Q_PROPERTY(bool earthJ3 MEMBER m_earthJ3 NOTIFY changed)
	Q_PROPERTY(bool earthJ4 MEMBER m_earthJ4 NOTIFY changed)
	Q_PROPERTY(bool rotatingEarthPole MEMBER m_rotatingEarthPole NOTIFY changed)


   public:
	using QObject::QObject;

	[[nodiscard]] bool enabled() const noexcept { return m_enabled; }

	[[nodiscard]] uint planetMask() const noexcept { return m_planetMask; }

	[[nodiscard]] bool useSb441() const noexcept { return m_useSb441; }

	[[nodiscard]] double integratorInitialStepDays() const noexcept { return m_integratorInitialStepDays; }

	[[nodiscard]] orbital_core::ForceModelSettings forceModel() const noexcept
	{
		return {static_cast<orbital_core::RelativityMode>(m_relativityMode), m_earthJ2, m_earthJ3, m_earthJ4,
		        m_rotatingEarthPole};
	}

	void setPlanetMask(uint v)
	{
		m_planetMask = v;
		emit changed();
	}

	void setUseSb441(bool v)
	{
		m_useSb441 = v;
		emit changed();
	}

	void save(QSettings& s) const;
	void load(QSettings& s);

	static constexpr uint kDefaultPlanetMask = static_cast<uint>(
	    (PLANET_BODY_MASK | MOON_BODY_MASK) & ~orbital_core::planetaryBodyMask(EphemerisPlanet::Pluto));

   signals:
	void changed();

   private:
	bool m_enabled{false};
	uint m_planetMask{kDefaultPlanetMask};
	bool m_useSb441{false};
	double m_integratorInitialStepDays{10.0};
	int m_relativityMode{0};
	bool m_earthJ2{true};
	bool m_earthJ3{true};
	bool m_earthJ4{true};
	bool m_rotatingEarthPole{true};
};
