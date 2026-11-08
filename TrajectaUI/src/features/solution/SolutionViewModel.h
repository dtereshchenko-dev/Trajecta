#pragma once

#include <QAbstractTableModel>
#include <QQmlEngine>
#include <QVariantMap>
#include <optional>
#include <orbital_core/Photometry.h>
#include <orbital_core/hazard/Moid.h>
#include <orbital_core/hazard/PhaClassifier.h>
#include "features/observations/ObservationTableModel.h"
#include "models/OrbitModel.h"
#include "models/DiagnosticsModel.h"
#include "CovarianceTableModel.h"
#include "features/orbit/OrbitSceneModel.h"

class SolutionViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(OrbitModel* orbit READ orbit CONSTANT)
	Q_PROPERTY(OrbitSceneModel* scene READ scene CONSTANT)
	Q_PROPERTY(ObservationTableModel* observations READ observations CONSTANT)
	Q_PROPERTY(CovarianceTableModel* covariance READ covariance CONSTANT)
	Q_PROPERTY(DiagnosticsModel* diagnostics READ diagnostics CONSTANT)
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(bool hasGaiaObservations READ hasGaiaObservations NOTIFY changed)
	Q_PROPERTY(bool hasGaiaFprObservations READ hasGaiaFprObservations NOTIFY changed)
	Q_PROPERTY(QString gaiaSummary READ gaiaSummary NOTIFY changed)
	Q_PROPERTY(bool converged READ converged NOTIFY changed)
	Q_PROPERTY(QString objectTitle READ objectTitle NOTIFY changed)
	Q_PROPERTY(QString designation READ designation NOTIFY changed)
	Q_PROPERTY(QString directory READ directory NOTIFY changed)
	Q_PROPERTY(QString epoch READ epoch NOTIFY changed)
	Q_PROPERTY(QString status READ status NOTIFY changed)
	Q_PROPERTY(int total READ total NOTIFY changed)
	Q_PROPERTY(int included READ included NOTIFY changed)
	Q_PROPERTY(double rms READ rms NOTIFY changed)
	Q_PROPERTY(double normalizedRms READ normalizedRms NOTIFY changed)
	Q_PROPERTY(double arcDays READ arcDays NOTIFY changed)
	Q_PROPERTY(bool cometary READ cometary NOTIFY changed)
	Q_PROPERTY(bool hasSigmas READ hasSigmas NOTIFY changed)
	Q_PROPERTY(QString epochMjd READ epochMjd NOTIFY changed)
	Q_PROPERTY(QString fitEpochUtc READ fitEpochUtc NOTIFY changed)
	Q_PROPERTY(QString calculationStartUtc READ calculationStartUtc NOTIFY changed)
	Q_PROPERTY(QString arcRange READ arcRange NOTIFY changed)
	Q_PROPERTY(int stationCount READ stationCount NOTIFY changed)
	Q_PROPERTY(QVariantList stateVector READ stateVector NOTIFY changed)
	Q_PROPERTY(double absoluteMagnitude READ absoluteMagnitude NOTIFY changed)
	Q_PROPERTY(double asteroidSlopeG READ asteroidSlopeG NOTIFY changed)
	Q_PROPERTY(bool cometPhotometry READ cometPhotometry NOTIFY changed)
	Q_PROPERTY(double cometTotalMagnitude READ cometTotalMagnitude NOTIFY changed)
	Q_PROPERTY(double cometNuclearMagnitude READ cometNuclearMagnitude NOTIFY changed)
	Q_PROPERTY(double cometK1 READ cometK1 NOTIFY changed)
	Q_PROPERTY(double cometK2 READ cometK2 NOTIFY changed)
	Q_PROPERTY(double cometNuclearPhaseCoefficient READ cometNuclearPhaseCoefficient NOTIFY changed)
	Q_PROPERTY(int magnitudeCount READ magnitudeCount NOTIFY changed)
	Q_PROPERTY(QVariantList sizeEstimates READ sizeEstimates NOTIFY changed)
	Q_PROPERTY(QVariantMap typicalSize READ typicalSize NOTIFY changed)
	Q_PROPERTY(double earthMoidKm READ earthMoidKm NOTIFY changed)
	Q_PROPERTY(QString phaStatus READ phaStatus NOTIFY changed)
	Q_PROPERTY(bool potentiallyHazardous READ potentiallyHazardous NOTIFY changed)
	Q_PROPERTY(QVariantList planetMoids READ planetMoids NOTIFY changed)
	Q_PROPERTY(QVariantList equinoctialElements READ equinoctialElements NOTIFY changed)
	Q_PROPERTY(bool nonGravEnabled READ nonGravEnabled NOTIFY changed)
	Q_PROPERTY(QString nonGravModel READ nonGravModel NOTIFY changed)
	Q_PROPERTY(QVariantList nonGravParameters READ nonGravParameters NOTIFY changed)

   public:
	explicit SolutionViewModel(QObject* parent = nullptr);
	bool publish(const orbital_core::FitPipelineResult& result, const QString& directory,
	             const orbital_core::hazard::PlanetMoidTable& moids = {});

	OrbitModel* orbit() { return &m_orbit; }

	OrbitSceneModel* scene() { return &m_scene; }

	ObservationTableModel* observations() { return &m_observations; }

	CovarianceTableModel* covariance() { return &m_covariance; }

	DiagnosticsModel* diagnostics() { return &m_diagnostics; }

	bool available() const { return m_snapshot.has_value(); }

	bool hasGaiaObservations() const;
	bool hasGaiaFprObservations() const;
	QString gaiaSummary() const;

	const std::optional<orbital_core::FitPipelineResult>& snapshot() const { return m_snapshot; }

	bool converged() const { return m_snapshot && m_snapshot->fit.converged; }

	QString objectTitle() const;
	QString designation() const;

	QString directory() const { return m_directory; }

	QString epoch() const;
	QString status() const;
	int total() const;
	int included() const;
	double rms() const;
	double normalizedRms() const;

	double arcDays() const { return m_arcDays; }

	bool cometary() const;
	bool hasSigmas() const;
	QString epochMjd() const;
	QString fitEpochUtc() const;
	QString calculationStartUtc() const;
	QString arcRange() const;
	int stationCount() const;
	QVariantList stateVector() const;
	double absoluteMagnitude() const;

	double asteroidSlopeG() const { return m_photometry.slopeG; }

	bool cometPhotometry() const;
	double cometTotalMagnitude() const;
	double cometNuclearMagnitude() const;

	double cometK1() const { return m_cometPhotometry.total.slope; }

	double cometK2() const { return m_cometPhotometry.nuclear.slope; }

	double cometNuclearPhaseCoefficient() const { return m_cometPhotometry.nuclear.phaseCoefficient; }

	const orbital_core::photometry::CometSolution& cometPhotometrySolution() const { return m_cometPhotometry; }

	int magnitudeCount() const { return static_cast<int>(m_photometry.count); }

	QVariantList sizeEstimates() const;
	QVariantMap typicalSize() const;
	double earthMoidKm() const;
	QString phaStatus() const;
	bool potentiallyHazardous() const;
	QVariantList planetMoids() const;

	QVariantList equinoctialElements() const { return m_equinoctialElements; }

	bool nonGravEnabled() const;
	QString nonGravModel() const;
	QVariantList nonGravParameters() const;

   signals:
	void changed();

   private:
	orbital_core::hazard::PhaStatus phaClassification() const;
	OrbitModel m_orbit;
	OrbitSceneModel m_scene;
	ObservationTableModel m_observations;
	CovarianceTableModel m_covariance;
	DiagnosticsModel m_diagnostics;
	std::optional<orbital_core::FitPipelineResult> m_snapshot;
	QString m_directory;
	double m_arcDays{};
	orbital_core::photometry::Solution m_photometry;
	orbital_core::photometry::CometSolution m_cometPhotometry;
	orbital_core::hazard::PlanetMoidTable m_moids;
	QVariantList m_equinoctialElements;
};
