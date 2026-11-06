#include "NonGravSettings.h"

void NonGravSettings::save(QSettings& s) const
{
	s.setValue("nonGravEnabled", m_enabled);
	s.setValue("nonGravModel", m_model);
	s.setValue("nonGravSolveMask", m_solveMask);
}

void NonGravSettings::load(QSettings& s)
{
	m_enabled = s.value("nonGravEnabled", false).toBool();
	m_model = s.value("nonGravModel", 1).toInt();
	m_solveMask = s.value("nonGravSolveMask", 0).toInt();
	emit changed();
}
