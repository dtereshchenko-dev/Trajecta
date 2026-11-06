#pragma once

#include <QObject>
#include <QList>
#include <QUrl>
#include <optional>
#include <QQmlEngine>
#include "controllers/OrbitEngine.h"
#include "features/solution/SolutionViewModel.h"
#include "features/ephemerides/EphemerisViewModel.h"
#include "features/encounters/CloseApproachesViewModel.h"
#include "features/risk/RiskMonitorViewModel.h"
#include "features/risk/MovViewModel.h"
#include "features/atmosphere/AtmosphereViewModel.h"
#include "features/settings/CalculationSettingsModel.h"
#include "features/settings/DataSourcesModel.h"
#include "features/observations/ObservationReviewModel.h"
#include "io/RecentFiles.h"

class ApplicationSession : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Created by the application")
	Q_PROPERTY(SolutionViewModel* solution READ solution CONSTANT)
	Q_PROPERTY(EphemerisViewModel* ephemerides READ ephemerides CONSTANT)
	Q_PROPERTY(CloseApproachesViewModel* approaches READ approaches CONSTANT)
	Q_PROPERTY(RiskMonitorViewModel* risk READ risk CONSTANT)
	Q_PROPERTY(MovViewModel* mov READ mov CONSTANT)
	Q_PROPERTY(AtmosphereViewModel* atmosphere READ atmosphere CONSTANT)
	Q_PROPERTY(RecentFiles* recentFiles READ recentFiles CONSTANT)
	Q_PROPERTY(CalculationSettingsModel* calculationSettings READ calculationSettings CONSTANT)
	Q_PROPERTY(DataSourcesModel* dataSources READ dataSources CONSTANT)
	Q_PROPERTY(ObservationReviewModel* observationReview READ observationReview CONSTANT)
	Q_PROPERTY(QString inputFile READ inputFile NOTIFY changed)
	Q_PROPERTY(QString inputName READ inputName NOTIFY changed)
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(bool cancelling READ cancelling NOTIFY changed)
	Q_PROPERTY(bool canRun READ canRun NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QString progress READ progress NOTIFY changed)
	Q_PROPERTY(QString theme READ theme WRITE setTheme NOTIFY themeChanged)

   public:
	ApplicationSession(FitSettingsAggregator& settings, QObject* parent = nullptr);

	SolutionViewModel* solution() { return &m_solution; }

	EphemerisViewModel* ephemerides() { return &m_ephemerides; }

	CloseApproachesViewModel* approaches() { return &m_approaches; }

	RiskMonitorViewModel* risk() { return &m_risk; }

	MovViewModel* mov() { return &m_mov; }

	AtmosphereViewModel* atmosphere() { return &m_atmosphere; }

	RecentFiles* recentFiles() { return &m_recentFiles; }

	CalculationSettingsModel* calculationSettings() { return &m_calculationSettings; }

	DataSourcesModel* dataSources() { return &m_dataSources; }

	ObservationReviewModel* observationReview() { return &m_observationReview; }

	QString inputFile() const { return m_inputFile; }

	QString inputName() const;

	bool busy() const
	{
		return m_busy || m_risk.busy() || m_ephemerides.busy() || m_atmosphere.busy() || m_approaches.busy() ||
		       m_mov.busy();
	}

	bool cancelling() const
	{
		return m_cancelling || m_risk.cancelling() || m_atmosphere.cancelling() || m_approaches.cancelling() ||
		       m_mov.cancelling();
	}

	bool canRun() const { return !busy() && !m_inputFile.isEmpty(); }

	QString error() const { return m_error; }

	QString progress() const
	{
		return m_mov.busy()           ? m_mov.progress()
		       : m_approaches.busy()  ? m_approaches.progress()
		       : m_atmosphere.busy()  ? m_atmosphere.progress()
		       : m_risk.busy()        ? m_risk.progress()
		       : m_ephemerides.busy() ? tr("Computing ephemeris...")
		                              : m_progress;
	}

	QString theme() const { return m_theme; }

	void setTheme(const QString& value);
	Q_INVOKABLE bool openInput(const QUrl& url);
	Q_INVOKABLE bool openRecent(int index);
	Q_INVOKABLE bool canDropInput(const QList<QUrl>& urls) const;
	Q_INVOKABLE void run();
	Q_INVOKABLE void cancel();
	Q_INVOKABLE void openDirectory();
	Q_INVOKABLE void dismissError();
	Q_INVOKABLE void acceptObservationReview();
	Q_INVOKABLE void cancelObservationReview();
	Q_INVOKABLE bool applySettings(bool recalculate = false);
	Q_INVOKABLE bool downloadGaiaFpr();

   signals:
	void changed();
	void themeChanged();
	void runStarted();
	void observationReviewReady();

   private:
	SolutionViewModel m_solution;
	EphemerisViewModel m_ephemerides;
	CloseApproachesViewModel m_approaches{m_solution};
	RiskMonitorViewModel m_risk{m_solution};
	AtmosphereViewModel m_atmosphere{m_solution, m_risk};
	RecentFiles m_recentFiles;
	QString m_inputFile;
	QString m_acceptedInputFile;
	bool m_loadRequested{true};
	QString m_error;
	QString m_progress;
	QString m_theme;
	bool m_busy{};
	bool m_cancelling{};
	OrbitEngine m_engine;
	MovViewModel m_mov{m_solution, [this] { return m_engine.currentSettings(); }};
	CalculationSettingsModel m_calculationSettings;
	DataSourcesModel m_dataSources;
	ObservationReviewModel m_observationReview;
	bool m_gaiaDownloading{};
	std::shared_ptr<orbital_core::ObservationSet> m_preparedInput;
	bool m_reviewingInput{};
};
