#include "FitSettingsAggregator.h"

#include <QSettings>

FitSettingsAggregator::FitSettingsAggregator(QObject* parent)
    : QObject(parent),
      m_perturbations(new PerturbationsSettings(this)),
      m_nonGrav(new NonGravSettings(this)),
      m_weighting(new WeightingSettings(this)),
      m_sigmaClip(new SigmaClipSettings(this)),
      m_photometry(new PhotometrySettings(this))
{
	connect(m_perturbations, &PerturbationsSettings::changed, this, &FitSettingsAggregator::anyChanged);
	connect(m_nonGrav, &NonGravSettings::changed, this, &FitSettingsAggregator::anyChanged);
	connect(m_weighting, &WeightingSettings::changed, this, &FitSettingsAggregator::anyChanged);
	connect(m_sigmaClip, &SigmaClipSettings::changed, this, &FitSettingsAggregator::anyChanged);
	connect(m_photometry, &PhotometrySettings::changed, this, &FitSettingsAggregator::anyChanged);
}

void FitSettingsAggregator::save()
{
	QSettings s;
	s.beginGroup("optimizer");
	m_perturbations->save(s);
	m_nonGrav->save(s);
	m_weighting->save(s);
	m_sigmaClip->save(s);
	m_photometry->save(s);
	s.endGroup();
}

void FitSettingsAggregator::load()
{
	QSettings s;
	s.beginGroup("optimizer");
	m_perturbations->load(s);
	m_nonGrav->load(s);
	m_weighting->load(s);
	m_sigmaClip->load(s);
	m_photometry->load(s);
	s.endGroup();
}
