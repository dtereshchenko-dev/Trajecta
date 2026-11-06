#include "WeightingSettings.h"

namespace {
	[[nodiscard]] int normalizedWeightProfile(int value) noexcept { return (value >= 1 && value <= 2) ? value : 1; }
}  // namespace

void WeightingSettings::save(QSettings& s) const
{
	s.setValue("useWeights", m_useWeights);
	s.setValue("observerSigmas", m_observerSigmas);
	s.setValue("timingEnabled", m_timingEnabled);
	s.setValue("debiasEnabled", m_debiasEnabled);
	s.setValue("weightProfile", normalizedWeightProfile(m_weightProfile));
}

void WeightingSettings::load(QSettings& s)
{
	m_useWeights = s.value("useWeights", true).toBool();
	m_observerSigmas = s.value("observerSigmas", true).toBool();
	m_timingEnabled = s.value("timingEnabled", false).toBool();
	m_debiasEnabled = s.value("debiasEnabled", true).toBool();
	m_weightProfile = normalizedWeightProfile(s.value("weightProfile", 1).toInt());
	emit changed();
}
