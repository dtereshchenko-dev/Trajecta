#pragma once

#include <QObject>
#include <QQmlEngine>
#include "settings/FitSettingsAggregator.h"

class CalculationSettingsModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(PerturbationsSettings* perturbations READ perturbations CONSTANT)
	Q_PROPERTY(WeightingSettings* weighting READ weighting CONSTANT)
	Q_PROPERTY(SigmaClipSettings* sigmaClip READ sigmaClip CONSTANT)
	Q_PROPERTY(NonGravSettings* nonGrav READ nonGrav CONSTANT)
	Q_PROPERTY(PhotometrySettings* photometry READ photometry CONSTANT)
	Q_PROPERTY(QString validationError READ validationError NOTIFY changed)
	Q_PROPERTY(bool valid READ valid NOTIFY changed)
	Q_PROPERTY(QString applicationVersion READ applicationVersion CONSTANT)
	Q_PROPERTY(QString qtVersion READ qtVersion CONSTANT)
	Q_PROPERTY(int nonGravDefaultMask READ nonGravDefaultMask NOTIFY changed)
	Q_PROPERTY(int nonGravSupportedMask READ nonGravSupportedMask NOTIFY changed)

   public:
	explicit CalculationSettingsModel(FitSettingsAggregator& settings, QObject* parent = nullptr);

	PerturbationsSettings* perturbations() { return &m_perturbations; }

	WeightingSettings* weighting() { return &m_weighting; }

	SigmaClipSettings* sigmaClip() { return &m_sigmaClip; }

	NonGravSettings* nonGrav() { return &m_nonGrav; }

	PhotometrySettings* photometry() { return &m_photometry; }

	QString validationError() const;

	bool valid() const { return validationError().isEmpty(); }

	QString applicationVersion() const;
	QString qtVersion() const;
	int nonGravDefaultMask() const;
	int nonGravSupportedMask() const;
	Q_INVOKABLE void reset();
	bool apply();

   signals:
	void changed();

   private:
	FitSettingsAggregator& m_settings;
	PerturbationsSettings m_perturbations;
	WeightingSettings m_weighting;
	SigmaClipSettings m_sigmaClip;
	NonGravSettings m_nonGrav;
	PhotometrySettings m_photometry;
};
