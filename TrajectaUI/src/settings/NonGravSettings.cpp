#include "NonGravSettings.h"

void NonGravSettings::save(QSettings& s) const
{
	s.setValue("enabled", m_enabled);
	s.setValue("model", m_model);
	s.setValue("solveMask", m_solveMask);
}

void NonGravSettings::load(QSettings& s)
{
	m_enabled = s.value("enabled", false).toBool();
	m_model = s.value("model", 1).toInt();
	m_solveMask = s.value("solveMask", 0).toInt();
	emit changed();
}
