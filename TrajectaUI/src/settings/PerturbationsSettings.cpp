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
	s.setValue("integrator/nBodyPropagation", m_enabled);
	s.setValue("integrator/initialStepDays", m_integratorInitialStepDays);
	s.setValue("perturbations/planetMask", m_planetMask);
	s.setValue("perturbations/useSb441", m_useSb441);
	s.setValue("forces/relativityMode", m_relativityMode);
	s.setValue("forces/earthJ2", m_earthJ2);
	s.setValue("forces/earthJ3", m_earthJ3);
	s.setValue("forces/earthJ4", m_earthJ4);
	s.setValue("forces/rotatingEarthPole", m_rotatingEarthPole);
}

void PerturbationsSettings::load(QSettings& s)
{
	m_enabled = s.value("integrator/nBodyPropagation", true).toBool();
	m_planetMask = s.value("perturbations/planetMask", kDefaultPlanetMask).toUInt();
	m_useSb441 = s.value("perturbations/useSb441", false).toBool();
	m_integratorInitialStepDays =
	    saneNBodyStepDays(s.value("integrator/initialStepDays", kDefaultNBodyStepDays).toDouble());
	const int relativity = s.value("forces/relativityMode", 0).toInt();
	m_relativityMode = relativity >= 0 && relativity <= 2 ? relativity : 0;
	s.remove("forces/sunJ2");
	m_earthJ2 = s.value("forces/earthJ2", true).toBool();
	m_earthJ3 = s.value("forces/earthJ3", true).toBool();
	m_earthJ4 = s.value("forces/earthJ4", true).toBool();
	m_rotatingEarthPole = s.value("forces/rotatingEarthPole", true).toBool();
	emit changed();
}
