#pragma once

#include "features/solution/SolutionViewModel.h"
#include "RiskEventTableModel.h"

#include <QObject>
#include <QQmlEngine>
#include <QString>
#include <cstdint>
#include <optional>
#include <utility>
#include <thread>

// Editable client draft. Epochs are UTC; admission freezes the accepted source
// revision and converts bounds to model TDB. No synthetic encounter event.
struct RiskMonitorDraft {
	enum class Method { Covariance, Observations };
	enum class Target { Earth, Moon, Mars };
	enum class IntervalMode { EventWindow, DateRange };
	Method method{};
	Target target{};
	IntervalMode intervalMode{};
	double startMjdUtc{};
	double endMjdUtc{};
	double altitudeKm{};
	int samples{};
	std::uint64_t seed{};
	quint64 sourceRevision{};
};

// Owns the draft, asynchronous immutable jobs and per-source run history.
// Workers never access mutable QML state; accepted source revisions guard
// publication. Full clone records live in the archive, table summaries in UI.
class RiskMonitorViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(int method READ method WRITE setMethod NOTIFY changed)
	Q_PROPERTY(int target READ target WRITE setTarget NOTIFY changed)
	Q_PROPERTY(int intervalMode READ intervalMode WRITE setIntervalMode NOTIFY changed)
	Q_PROPERTY(QString centerUtc READ centerUtc WRITE setCenterUtc NOTIFY changed)
	Q_PROPERTY(double halfWindowHours READ halfWindowHours WRITE setHalfWindowHours NOTIFY changed)
	Q_PROPERTY(QString rangeStartUtc READ rangeStartUtc WRITE setRangeStartUtc NOTIFY changed)
	Q_PROPERTY(QString rangeEndUtc READ rangeEndUtc WRITE setRangeEndUtc NOTIFY changed)
	Q_PROPERTY(double rangeYears READ rangeYears WRITE setRangeYears NOTIFY changed)
	Q_PROPERTY(QString samples READ samples WRITE setSamples NOTIFY changed)
	Q_PROPERTY(QString seed READ seed WRITE setSeed NOTIFY changed)
	Q_PROPERTY(int boundary READ boundary WRITE setBoundary NOTIFY changed)
	Q_PROPERTY(double altitudeKm READ altitudeKm NOTIFY changed)
	Q_PROPERTY(QString sourceError READ sourceError NOTIFY changed)
	Q_PROPERTY(QString formError READ formError NOTIFY changed)
	Q_PROPERTY(bool inputReady READ inputReady NOTIFY changed)
	Q_PROPERTY(QString fitEpoch READ fitEpoch NOTIFY changed)
	Q_PROPERTY(QString covarianceSize READ covarianceSize NOTIFY changed)
	Q_PROPERTY(int observationCount READ observationCount NOTIFY changed)
	Q_PROPERTY(int radarCount READ radarCount NOTIFY changed)
	Q_PROPERTY(QString windowDuration READ windowDuration NOTIFY changed)
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(bool cancelling READ cancelling NOTIFY changed)
	Q_PROPERTY(bool canRun READ canRun NOTIFY changed)
	Q_PROPERTY(QString progress READ progress NOTIFY changed)
	Q_PROPERTY(double progressFraction READ progressFraction NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QStringList runs READ runs NOTIFY changed)
	Q_PROPERTY(int selectedRun READ selectedRun WRITE setSelectedRun NOTIFY changed)
	Q_PROPERTY(QString summary READ summary NOTIFY changed)
	Q_PROPERTY(QString runDetails READ runDetails NOTIFY changed)
	Q_PROPERTY(QString warning READ warning NOTIFY changed)
	Q_PROPERTY(QString archivePath READ archivePath NOTIFY changed)
	Q_PROPERTY(RiskEventTableModel* events READ events CONSTANT)

   public:
	explicit RiskMonitorViewModel(SolutionViewModel& solution, QObject* parent = nullptr);
	~RiskMonitorViewModel() override;

	bool available() const { return m_available; }

	int method() const { return m_method; }

	void setMethod(int value);

	int target() const { return m_target; }

	void setTarget(int value);

	int intervalMode() const { return m_intervalMode; }

	void setIntervalMode(int value);

	QString centerUtc() const { return m_centerUtc; }

	void setCenterUtc(const QString& value);

	double halfWindowHours() const { return m_halfWindowHours; }

	void setHalfWindowHours(double value);

	QString rangeStartUtc() const { return m_rangeStartUtc; }

	void setRangeStartUtc(const QString& value);

	QString rangeEndUtc() const { return m_rangeEndUtc; }

	void setRangeEndUtc(const QString& value);
	double rangeYears() const;
	void setRangeYears(double value);

	QString samples() const { return m_samples[m_method]; }

	void setSamples(const QString& value);

	QString seed() const { return m_seed; }

	void setSeed(const QString& value);

	int boundary() const { return m_target == 0 ? m_boundary : 0; }

	void setBoundary(int value);

	double altitudeKm() const { return boundary() == 1 ? 100.0 : 0.0; }

	QString sourceError() const;
	QString formError() const;

	bool inputReady() const { return sourceError().isEmpty() && formError().isEmpty(); }

	QString fitEpoch() const;
	QString covarianceSize() const;

	int observationCount() const { return m_observationCount; }

	int radarCount() const { return m_radarCount; }

	QString windowDuration() const;
	std::optional<RiskMonitorDraft> draft() const;

	bool busy() const { return m_busy; }

	bool cancelling() const { return m_cancelling; }

	bool canRun() const { return !m_busy && !m_blocked && inputReady(); }

	void setBlocked(bool value);

	QString progress() const { return m_progress; }

	double progressFraction() const { return m_progressFraction; }

	QString error() const { return m_error; }

	QStringList runs() const;
	QVariantList atmosphericSources() const;

	int selectedRun() const { return m_selectedRun; }

	void setSelectedRun(int value);
	QString summary() const;
	QString runDetails() const;
	QString warning() const;
	QString archivePath() const;

	RiskEventTableModel* events() { return &m_events; }

	Q_INVOKABLE void run();
	Q_INVOKABLE void cancel();
	Q_INVOKABLE void openArchive(bool report = false);

   signals:
	void changed();
	void runStarted();

   private:
	void refreshSource();
	std::optional<std::pair<double, double>> intervalBounds() const;
	SolutionViewModel& m_solution;
	bool m_available{};
	QString m_covarianceError;
	QString m_observationError;
	double m_epochMjdTdb{};
	int m_parameters{};
	int m_observationCount{};
	int m_radarCount{};
	quint64 m_sourceRevision{};
	QString m_sourceObject;
	int m_method{};
	int m_target{};
	int m_intervalMode{};
	int m_boundary{1};
	QString m_centerUtc;
	double m_halfWindowHours{24.0};
	bool m_centerInitialized{};
	QString m_rangeStartUtc;
	QString m_rangeEndUtc;
	bool m_rangeInitialized{};
	QString m_samples[2]{QStringLiteral("10000"), QStringLiteral("1000")};
	QString m_seed{QStringLiteral("42")};

	struct Run {
		std::shared_ptr<orbital_core::hazard::IntervalMcResult> result;
		QString label;
		QString archive;
		double seconds{};
	};

	RiskEventTableModel m_events;
	std::vector<Run> m_runs;
	int m_selectedRun{-1};
	bool m_busy{};
	bool m_cancelling{};
	bool m_blocked{};
	QString m_progress;
	QString m_error;
	double m_progressFraction{};
	std::jthread m_worker;
};
