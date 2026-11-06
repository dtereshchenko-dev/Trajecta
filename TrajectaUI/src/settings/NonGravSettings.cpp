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
	if (m_model == 5) {
		m_model = 3;
		if (m_solveMask == 0)
			m_solveMask = 15;
	}
	emit changed();
}
