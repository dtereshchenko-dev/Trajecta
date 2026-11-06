#pragma once
#include "features/solution/SolutionViewModel.h"
#include "controllers/OrbitWorker.h"
#include <functional>
#include <thread>

class MovViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(bool cancelling READ cancelling NOTIFY changed)
	Q_PROPERTY(bool canRun READ canRun NOTIFY changed)
	Q_PROPERTY(QString sourceName READ sourceName NOTIFY changed)
	Q_PROPERTY(QString seedUtc READ seedUtc NOTIFY changed)
	Q_PROPERTY(bool externalSource READ externalSource NOTIFY changed)
	Q_PROPERTY(int observationCount READ observationCount NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QString progress READ progress NOTIFY changed)
	Q_PROPERTY(QVariantList runs READ runs NOTIFY changed)
	Q_PROPERTY(int selectedRun READ selectedRun WRITE setSelectedRun NOTIFY changed)
	Q_PROPERTY(QVariantMap result READ result NOTIFY changed)
   public:
	MovViewModel(SolutionViewModel& solution, std::function<OptimizationSettings()> settings,
	             QObject* parent = nullptr);
	~MovViewModel() override;

	bool busy() const { return m_busy; }

	bool cancelling() const { return m_cancelling; }

	bool canRun() const { return !m_busy && !m_blocked && observationCount() >= 3; }

	bool externalSource() const { return m_external.has_value(); }

	QString sourceName() const;
	QString seedUtc() const;
	int observationCount() const;

	QString error() const { return m_error; }

	QString progress() const { return m_progress; }

	QVariantList runs() const { return m_runs; }

	int selectedRun() const { return m_selected; }

	QVariantMap result() const;
	void setSelectedRun(int index);
	void setBlocked(bool blocked);
	Q_INVOKABLE void selectFile(const QUrl& url);
	Q_INVOKABLE void useAcceptedSolution();
	Q_INVOKABLE void run(const QVariantMap& options);
	Q_INVOKABLE void cancel();
	Q_INVOKABLE void openReport();
	Q_INVOKABLE void openArchive();
   signals:
	void changed();
	void runStarted();
	void runFailed();

   private:
	SolutionViewModel& m_solution;
	std::function<OptimizationSettings()> m_settings;
	std::optional<orbital_core::ObservationSet> m_external;
	QString m_file, m_progress, m_error;
	QVariantList m_runs;
	int m_selected{-1};
	bool m_busy{}, m_blocked{}, m_cancelling{};
	quint64 m_revision{};
	std::jthread m_worker;
	void resetSource();
};
