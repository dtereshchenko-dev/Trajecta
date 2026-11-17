#include "CalculationSettingsModel.h"
#include "orbital_core/RuntimeData.h"

#include <QMetaProperty>
#include <QSignalBlocker>
#include <cmath>
#include <orbital_core/nongrav/NonGravModel.h>

namespace {
	template <typename T>
	void copySettings(const T& source, T& target)
	{
		const auto& meta = T::staticMetaObject;
		for (int i = meta.propertyOffset(); i < meta.propertyCount(); ++i) {
			const auto property = meta.property(i);
			if (property.isReadable() && property.isWritable())
				property.write(&target, property.read(&source));
		}
	}
}  // namespace

CalculationSettingsModel::CalculationSettingsModel(FitSettingsAggregator& settings, QObject* parent)
    : QObject(parent), m_settings(settings)
{
	connect(&m_perturbations, &PerturbationsSettings::changed, this, &CalculationSettingsModel::changed);
	connect(&m_weighting, &WeightingSettings::changed, this, &CalculationSettingsModel::changed);
	connect(&m_sigmaClip, &SigmaClipSettings::changed, this, &CalculationSettingsModel::changed);
	connect(&m_nonGrav, &NonGravSettings::changed, this, &CalculationSettingsModel::changed);
	connect(&m_photometry, &PhotometrySettings::changed, this, &CalculationSettingsModel::changed);
	reset();
}

void CalculationSettingsModel::reset()
{
	copySettings(*m_settings.perturbations(), m_perturbations);
	copySettings(*m_settings.weighting(), m_weighting);
	copySettings(*m_settings.sigmaClip(), m_sigmaClip);
	copySettings(*m_settings.nonGrav(), m_nonGrav);
	copySettings(*m_settings.photometry(), m_photometry);
	emit changed();
}

QString CalculationSettingsModel::validationError() const
{
	const auto step = m_perturbations.integratorInitialStepDays();
	const auto photometry = m_photometry.options();
	if (!std::isfinite(photometry.asteroidG) || !std::isfinite(photometry.cometK1) ||
	    !std::isfinite(photometry.cometK2) || !std::isfinite(photometry.cometNuclearPhaseCoefficient))
		return tr("Photometry: all coefficients must be finite numbers.");
	if (photometry.asteroidG < -99.999 || photometry.asteroidG > 999.999 ||
	    std::abs(photometry.asteroidG * 1000.0 - std::round(photometry.asteroidG * 1000.0)) > 1e-8)
		return tr("Photometry: G must fit the OEF format and have at most three decimal places.");
	if (!std::isfinite(step) || step < 1e-6)
		return tr("Integrator: initial step must be at least 0.000001 days.");
	if (m_perturbations.enabled() && !orbital_core::planetaryEphemerisAvailable())
		return tr("Integrator: DE440 is unavailable.");
	if (m_perturbations.enabled() && m_perturbations.useSb441() && !orbital_core::asteroidEphemerisAvailable())
		return tr("Perturbations: SB441 is unavailable.");
	const int relativity = static_cast<int>(m_perturbations.forceModel().relativity);
	if (relativity < 0 || relativity > 2)
		return tr("Forces: select a valid relativity model.");
	if (m_weighting.weightProfile() < 1 || m_weighting.weightProfile() > 2)
		return tr("Fit: select a valid weighting profile.");
	if (!std::isfinite(m_sigmaClip.threshold()) || m_sigmaClip.threshold() <= 0 ||
	    m_sigmaClip.threshold() >= 60.0 / 1.1)
		return tr("Fit: rejection threshold must be positive and below 54.545 sigma.");
	if (m_nonGrav.model() < 1 || m_nonGrav.model() > 4)
		return tr("Non-gravitational: select a valid model.");
	const int supportedMask = nonGravSupportedMask();
	if (m_nonGrav.solveMask() < 0 || (m_nonGrav.solveMask() & ~supportedMask) != 0)
		return tr("Non-gravitational: selected coefficients are not supported by this model.");
	if (m_nonGrav.solveMask() != 0 &&
	    m_nonGrav.solveMask() !=
	        static_cast<int>(NonGravModel::sanitizeFitParameterMaskFor(
	            static_cast<NgModel>(m_nonGrav.model()), static_cast<NgParameterMask>(m_nonGrav.solveMask()))))
		return tr("Non-gravitational: select a supported coefficient group.");
	return {};
}

bool CalculationSettingsModel::apply()
{
	if (!valid())
		return false;
	{
		const QSignalBlocker blocker(&m_settings);
		copySettings(m_perturbations, *m_settings.perturbations());
		copySettings(m_weighting, *m_settings.weighting());
		copySettings(m_sigmaClip, *m_settings.sigmaClip());
		copySettings(m_nonGrav, *m_settings.nonGrav());
		copySettings(m_photometry, *m_settings.photometry());
	}
	m_settings.save();
	emit m_settings.anyChanged();
	return true;
}

QString CalculationSettingsModel::applicationVersion() const { return QStringLiteral(TRAJECTA_VERSION); }

QString CalculationSettingsModel::qtVersion() const { return QString::fromLatin1(qVersion()); }

int CalculationSettingsModel::nonGravDefaultMask() const
{
	return static_cast<int>(NonGravModel::defaultFitParameterMaskFor(static_cast<NgModel>(m_nonGrav.model())));
}

int CalculationSettingsModel::nonGravSupportedMask() const
{
	return static_cast<int>(NonGravModel::sanitizeFitParameterMaskFor(static_cast<NgModel>(m_nonGrav.model()),
	                                                                  static_cast<NgParameterMask>(15)));
}
