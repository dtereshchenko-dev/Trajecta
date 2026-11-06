#pragma once

#include <QObject>
#include <QQmlEngine>

#include "NonGravSettings.h"
#include "PerturbationsSettings.h"
#include "SigmaClipSettings.h"
#include "WeightingSettings.h"
#include "PhotometrySettings.h"

class FitSettingsAggregator : public QObject {
	Q_OBJECT

	Q_PROPERTY(PerturbationsSettings* perturbations READ perturbations CONSTANT)
	Q_PROPERTY(NonGravSettings* nonGrav READ nonGrav CONSTANT)
	Q_PROPERTY(WeightingSettings* weighting READ weighting CONSTANT)
	Q_PROPERTY(SigmaClipSettings* sigmaClip READ sigmaClip CONSTANT)
	Q_PROPERTY(PhotometrySettings* photometry READ photometry CONSTANT)

   public:
	explicit FitSettingsAggregator(QObject* parent = nullptr);

	[[nodiscard]] PerturbationsSettings* perturbations() const noexcept { return m_perturbations; }

	[[nodiscard]] NonGravSettings* nonGrav() const noexcept { return m_nonGrav; }

	[[nodiscard]] WeightingSettings* weighting() const noexcept { return m_weighting; }

	[[nodiscard]] SigmaClipSettings* sigmaClip() const noexcept { return m_sigmaClip; }

	[[nodiscard]] PhotometrySettings* photometry() const noexcept { return m_photometry; }

	Q_INVOKABLE void save();
	Q_INVOKABLE void load();

   signals:
	void anyChanged();

   private:
	PerturbationsSettings* m_perturbations{nullptr};
	NonGravSettings* m_nonGrav{nullptr};
	WeightingSettings* m_weighting{nullptr};
	SigmaClipSettings* m_sigmaClip{nullptr};
	PhotometrySettings* m_photometry{nullptr};
};
