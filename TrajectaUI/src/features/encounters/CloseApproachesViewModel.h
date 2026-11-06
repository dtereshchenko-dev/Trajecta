#pragma once

#include "CloseApproachTableModel.h"
#include "features/solution/SolutionViewModel.h"
#include <thread>

class CloseApproachesViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(CloseApproachTableModel* table READ table CONSTANT)
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(bool cancelling READ cancelling NOTIFY changed)
	Q_PROPERTY(bool canRun READ canRun NOTIFY changed)
	Q_PROPERTY(bool hasResult READ hasResult NOTIFY changed)
	Q_PROPERTY(QString progress READ progress NOTIFY changed)
	Q_PROPERTY(QString summary READ summary NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QString archivePath READ archivePath NOTIFY changed)
   public:
	explicit CloseApproachesViewModel(SolutionViewModel& solution, QObject* parent = nullptr);
	~CloseApproachesViewModel() override;

	CloseApproachTableModel* table() { return &m_table; }

	bool available() const;

	bool busy() const { return m_busy; }

	bool cancelling() const { return m_cancelling; }

	bool canRun() const { return available() && !m_busy && !m_blocked; }

	bool hasResult() const { return m_hasResult; }

	QString progress() const { return m_progress; }

	QString summary() const { return m_summary; }

	QString error() const { return m_error; }

	QString archivePath() const { return m_archive; }

	void setBlocked(bool value);
	Q_INVOKABLE QString validationError(const QString& startUtc, const QString& endUtc, double maximumAu) const;
	Q_INVOKABLE void run(const QString& startUtc, const QString& endUtc, double maximumAu);
	Q_INVOKABLE void cancel();
	Q_INVOKABLE void openArchive();
   signals:
	void changed();
	void completed();

   private:
	SolutionViewModel& m_solution;
	CloseApproachTableModel m_table;
	std::jthread m_worker;
	quint64 m_revision{};
	bool m_busy{}, m_cancelling{}, m_blocked{}, m_hasResult{};
	QString m_progress, m_summary, m_error, m_archive;
};
