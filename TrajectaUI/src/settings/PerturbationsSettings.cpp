#include "PerturbationsSettings.h"

#include <cmath>

namespace {
	constexpr double kDefaultNBodyStepDays = 10.0;
	constexpr double kMinSaneNBodyStepDays = 1.0e-6;

	double saneNBodyStepDays(double stepDays)
	{
		if (!std::isfinite(stepDays) || stepDays < kMinSaneNBodyStepDays)
			return kDefaultNBodyStepDays;
		return stepDays;
	}
}  // namespace

void PerturbationsSettings::save(QSettings& s) const
{
	s.setValue("pertEnabled", m_enabled);
	s.setValue("planetMask", m_planetMask);
	s.setValue("useSb441", m_useSb441);
	s.setValue("integratorInitialStepDays", m_integratorInitialStepDays);
}

void PerturbationsSettings::load(QSettings& s)
{
	m_enabled = s.value("pertEnabled", true).toBool();
	m_planetMask = s.value("planetMask", kDefaultPlanetMask).toUInt();
	m_useSb441 = s.value("useSb441", false).toBool();
	m_integratorInitialStepDays =
	    saneNBodyStepDays(s.value("integratorInitialStepDays", kDefaultNBodyStepDays).toDouble());
	emit changed();
}
