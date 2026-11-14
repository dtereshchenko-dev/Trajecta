#include "FitSettingsAggregator.h"

#include <QMap>
#include <QSettings>

namespace {
	void migrateOptimizerSettings(QSettings& settings)
	{
		settings.beginGroup("optimizer");
		const QStringList keys = settings.allKeys();
		settings.endGroup();
		if (keys.isEmpty())
			return;

		const QMap<QString, QString> destinations{
		    {"pertEnabled", "integrator/nBodyPropagation"},
		    {"integratorInitialStepDays", "integrator/initialStepDays"},
		    {"planetMask", "perturbations/planetMask"},
		    {"useSb441", "perturbations/useSb441"},
		    {"nonGravEnabled", "non-gravitational/enabled"},
		    {"nonGravModel", "non-gravitational/model"},
		    {"nonGravSolveMask", "non-gravitational/solveMask"},
		    {"useWeights", "fit/useWeights"},
		    {"observerSigmas", "fit/observerSigmas"},
		    {"timingEnabled", "fit/timingEnabled"},
		    {"debiasEnabled", "fit/debiasEnabled"},
		    {"weightProfile", "fit/weightProfile"},
		    {"sigmaClipEnabled", "fit/sigmaClipEnabled"},
		    {"sigmaClipThreshold", "fit/sigmaClipThreshold"},
		    {"photometry/asteroidG", "photometry/asteroidG"},
		    {"photometry/cometK1", "photometry/cometK1"},
		    {"photometry/cometK2", "photometry/cometK2"},
		    {"photometry/nuclearPhaseCoefficient", "photometry/nuclearPhaseCoefficient"},
		};

		for (const QString& key : keys) {
			const auto destination = destinations.constFind(key);
			if (destination != destinations.cend() && !settings.contains(destination.value()))
				settings.setValue(destination.value(), settings.value("optimizer/" + key));
		}
		settings.sync();
		if (settings.status() != QSettings::NoError) {
			qWarning("Could not save migrated settings");
			return;
		}
		settings.remove("optimizer");
		settings.sync();
		if (settings.status() != QSettings::NoError)
			qWarning("Could not remove old settings section");
	}
}  // namespace

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
	m_perturbations->save(s);
	s.beginGroup("non-gravitational");
	m_nonGrav->save(s);
	s.endGroup();
	s.beginGroup("fit");
	m_weighting->save(s);
	m_sigmaClip->save(s);
	s.endGroup();
	m_photometry->save(s);
}

void FitSettingsAggregator::load()
{
	QSettings s;
	s.setFallbacksEnabled(false);
	migrateOptimizerSettings(s);
	m_perturbations->load(s);
	s.beginGroup("non-gravitational");
	m_nonGrav->load(s);
	s.endGroup();
	s.beginGroup("fit");
	m_weighting->load(s);
	m_sigmaClip->load(s);
	s.endGroup();
	m_photometry->load(s);
}
