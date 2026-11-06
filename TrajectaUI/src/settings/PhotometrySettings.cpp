#include "PhotometrySettings.h"

void PhotometrySettings::save(QSettings& settings) const
{
	settings.setValue("photometry/asteroidG", m_asteroidG);
	settings.setValue("photometry/cometK1", m_cometK1);
	settings.setValue("photometry/cometK2", m_cometK2);
	settings.setValue("photometry/nuclearPhaseCoefficient", m_nuclearPhaseCoefficient);
}

void PhotometrySettings::load(QSettings& settings)
{
	const orbital_core::PhotometryOptions defaults;
	m_asteroidG = settings.value("photometry/asteroidG", defaults.asteroidG).toDouble();
	m_cometK1 = settings.value("photometry/cometK1", defaults.cometK1).toDouble();
	m_cometK2 = settings.value("photometry/cometK2", defaults.cometK2).toDouble();
	m_nuclearPhaseCoefficient =
	    settings.value("photometry/nuclearPhaseCoefficient", defaults.cometNuclearPhaseCoefficient).toDouble();
	emit changed();
}
