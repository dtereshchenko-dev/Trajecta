#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QVariantMap>
#include <QUrl>
#include "features/solution/SolutionViewModel.h"
#include "features/risk/RiskMonitorViewModel.h"
#include "strewn_core/AtmosphericInput.h"
#include "AtmosphereArchive.h"
#include <thread>

struct AtmosphereParameters {
	strewn::PhysicsConfig physics;
	strewn::IntegrationOptions integration;
	strewn::MonteCarloConfig monteCarlo;
	strewn::MapOptions maps;
	std::variant<strewn::AutomaticEnvironment, strewn::EnvironmentFile, strewn::IdealizedEnvironment> environment;
};

class AtmosphereViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(QVariantMap values READ values NOTIFY changed)
	Q_PROPERTY(int sourceMode READ sourceMode WRITE setSourceMode NOTIFY changed)
	Q_PROPERTY(int sizeMode READ sizeMode WRITE setSizeMode NOTIFY changed)
	Q_PROPERTY(int environmentMode READ environmentMode WRITE setEnvironmentMode NOTIFY changed)
	Q_PROPERTY(int preset READ preset NOTIFY changed)
	Q_PROPERTY(bool offline READ offline WRITE setOffline NOTIFY changed)
	Q_PROPERTY(bool convergence READ convergence WRITE setConvergence NOTIFY changed)
	Q_PROPERTY(int mapMethod READ mapMethod WRITE setMapMethod NOTIFY changed)
	Q_PROPERTY(QString archivePath READ archivePath NOTIFY changed)
	Q_PROPERTY(QString environmentPath READ environmentPath NOTIFY changed)
	Q_PROPERTY(QVariantList sources READ sources NOTIFY changed)
	Q_PROPERTY(QString formError READ formError NOTIFY changed)
	Q_PROPERTY(QString sourceError READ sourceError NOTIFY changed)
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(bool loading READ loading NOTIFY changed)
	Q_PROPERTY(bool cancelling READ cancelling NOTIFY changed)
	Q_PROPERTY(bool canRun READ canRun NOTIFY changed)
	Q_PROPERTY(QString progress READ progress NOTIFY changed)
	Q_PROPERTY(double progressFraction READ progressFraction NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QStringList events READ events NOTIFY changed)
	Q_PROPERTY(int selectedEvent READ selectedEvent WRITE setSelectedEvent NOTIFY changed)
	Q_PROPERTY(QString sourceSummary READ sourceSummary NOTIFY changed)
	Q_PROPERTY(QStringList runs READ runs NOTIFY changed)
	Q_PROPERTY(int selectedRun READ selectedRun WRITE setSelectedRun NOTIFY changed)
	Q_PROPERTY(QString summary READ summary NOTIFY changed)
	Q_PROPERTY(QString runError READ runError NOTIFY changed)
	Q_PROPERTY(QString resultPath READ resultPath NOTIFY changed)
	Q_PROPERTY(bool reportAvailable READ reportAvailable NOTIFY changed)

   public:
	AtmosphereViewModel(SolutionViewModel& solution, RiskMonitorViewModel& risk, QObject* parent = nullptr);
	~AtmosphereViewModel() override;

	QVariantMap values() const { return m_values; }

	int sourceMode() const { return m_sourceMode; }

	void setSourceMode(int value);

	int sizeMode() const { return m_sizeMode; }

	void setSizeMode(int value);

	int environmentMode() const { return m_environmentMode; }

	void setEnvironmentMode(int value);
	int preset() const;

	bool offline() const { return m_offline; }

	void setOffline(bool value);

	bool convergence() const { return m_convergence; }

	void setConvergence(bool value);

	int mapMethod() const { return m_mapMethod; }

	void setMapMethod(int value);

	QString archivePath() const { return m_archivePath; }

	QString environmentPath() const { return m_environmentPath; }

	QVariantList sources() const { return m_sources; }

	QString formError() const;
	QString sourceError() const;

	bool busy() const { return m_busy || m_loading; }

	bool loading() const { return m_loading; }

	bool cancelling() const { return m_cancelling; }

	bool canRun() const { return !busy() && !m_blocked && sourceError().isEmpty() && formError().isEmpty(); }

	void setBlocked(bool value);

	QString progress() const { return m_progress; }

	double progressFraction() const { return m_progressFraction; }

	QString error() const { return m_error; }

	QStringList events() const;

	int selectedEvent() const { return m_selectedEvent; }

	void setSelectedEvent(int value);
	QString sourceSummary() const;
	QStringList runs() const;

	int selectedRun() const { return m_selectedRun; }

	void setSelectedRun(int value);
	QString summary() const;
	QString runError() const;
	QString resultPath() const;
	bool reportAvailable() const;
	AtmosphereParameters parameterOptions() const;
	Q_INVOKABLE void setValue(const QString& key, const QString& value);
	Q_INVOKABLE void selectPreset(int index);
	Q_INVOKABLE void selectArchive(const QUrl& url);
	Q_INVOKABLE void selectEnvironment(const QUrl& url);
	Q_INVOKABLE void run();
	Q_INVOKABLE void cancel();
	Q_INVOKABLE void openResult(bool report = false);

   signals:
	void changed();
	void runStarted();

   private:
	void readMaterial(const strewn::JobSpec& job);
	void refreshMagnitude();
	void refreshSource();
	SolutionViewModel& m_solution;
	RiskMonitorViewModel& m_risk;
	QVariantMap m_values;
	QVariantList m_sources;
	QString m_sourceObject;
	quint64 m_sourceRevision{};
	QString m_archivePath;
	QString m_environmentPath;
	bool m_magnitudeEdited{};
	bool m_offline{};
	bool m_convergence{true};
	int m_sourceMode{};
	int m_sizeMode{};
	int m_environmentMode{};
	int m_mapMethod{};

	struct Run {
		QString label;
		QString path;
		QString summary;
		QString error;
	};

	std::shared_ptr<const atmosphere::EntryArchive> m_archive;
	QString m_archiveError;
	int m_selectedEvent{-1};
	std::vector<Run> m_runs;
	int m_selectedRun{-1};
	bool m_busy{};
	bool m_loading{};
	bool m_cancelling{};
	bool m_blocked{};
	QString m_progress;
	QString m_error;
	double m_progressFraction{};
	std::jthread m_worker;
};
