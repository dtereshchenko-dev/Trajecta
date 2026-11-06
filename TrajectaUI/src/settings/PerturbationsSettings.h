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


   public:
	using QObject::QObject;

	[[nodiscard]] bool enabled() const noexcept { return m_enabled; }

	[[nodiscard]] uint planetMask() const noexcept { return m_planetMask; }

	[[nodiscard]] bool useSb441() const noexcept { return m_useSb441; }

	[[nodiscard]] double integratorInitialStepDays() const noexcept { return m_integratorInitialStepDays; }

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
};
