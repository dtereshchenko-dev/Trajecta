#include "SigmaClipSettings.h"

void SigmaClipSettings::save(QSettings& s) const
{
	s.setValue("sigmaClipEnabled", m_enabled);
	s.setValue("sigmaClipThreshold", m_threshold);
}

void SigmaClipSettings::load(QSettings& s)
{
	m_enabled = s.value("sigmaClipEnabled", true).toBool();
	m_threshold = s.value("sigmaClipThreshold", 3.0).toDouble();
	emit changed();
}
