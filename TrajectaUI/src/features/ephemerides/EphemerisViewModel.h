#pragma once

#include "EphemerisTableModel.h"

#include <QObject>
#include <QQmlEngine>
#include <QString>
#include <QStringList>
#include <QUrl>
#include <orbital_core/Ephemeris.h>
#include <orbital_core/FitPipeline.h>
#include <orbital_core/Photometry.h>

#include <memory>
#include <optional>
#include <thread>

// Ephemerides page of the workspace. Holds the accepted solution (fit-epoch
// state with the effective force model of that fit, non-gravitational
// parameters included), the form of the page and the last computed table.
// orbital_core::computeEphemeris runs on a worker thread; a newer request or
// cancel() stops the previous one and its late result is discarded. The form
// is interpreted in the selected time scale: both dates are combined with the
// start time from the page settings and the end epoch is inclusive. Exports:
// JSON carries every numeric field of a row at full precision, TXT the
// table as displayed with a description of the request.
class EphemerisViewModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(EphemerisTableModel* table READ table CONSTANT)
	Q_PROPERTY(bool available READ available NOTIFY changed)
	Q_PROPERTY(bool busy READ busy NOTIFY changed)
	Q_PROPERTY(QString error READ error NOTIFY changed)
	Q_PROPERTY(QString summary READ summary NOTIFY changed)
	Q_PROPERTY(int epochCount READ epochCount NOTIFY changed)
	Q_PROPERTY(int plannedEpochs READ plannedEpochs NOTIFY formChanged)
	Q_PROPERTY(QString formError READ formError NOTIFY formChanged)
	Q_PROPERTY(bool startValid READ startValid NOTIFY formChanged)
	Q_PROPERTY(bool endValid READ endValid NOTIFY formChanged)
	Q_PROPERTY(bool timeValid READ timeValid NOTIFY formChanged)
	Q_PROPERTY(bool canGenerate READ canGenerate NOTIFY changed)
	Q_PROPERTY(int layout READ layout WRITE setLayout NOTIFY formChanged)
	Q_PROPERTY(int center READ center WRITE setCenter NOTIFY formChanged)
	Q_PROPERTY(QString station READ station WRITE setStation NOTIFY formChanged)
	Q_PROPERTY(QString startDate READ startDate WRITE setStartDate NOTIFY formChanged)
	Q_PROPERTY(QString endDate READ endDate WRITE setEndDate NOTIFY formChanged)
	Q_PROPERTY(QString startTime READ startTime WRITE setStartTime NOTIFY formChanged)
	Q_PROPERTY(double stepDays READ stepDays WRITE setStepDays NOTIFY formChanged)
	Q_PROPERTY(int plane READ plane WRITE setPlane NOTIFY formChanged)
	Q_PROPERTY(int scale READ scale WRITE setScale NOTIFY formChanged)
	Q_PROPERTY(bool lightTime READ lightTime WRITE setLightTime NOTIFY formChanged)
	Q_PROPERTY(QString scaleName READ scaleName NOTIFY formChanged)
	Q_PROPERTY(QStringList exportFormats READ exportFormats CONSTANT)

   public:
	enum Center { Heliocentric = 0, Geocentric = 1, Topocentric = 2 };
	Q_ENUM(Center)

	explicit EphemerisViewModel(QObject* parent = nullptr);
	~EphemerisViewModel() override;

	void setSolution(const orbital_core::FitPipelineResult& result, const QString& designation,
	                 std::optional<double> absoluteMagnitude,
	                 const orbital_core::photometry::CometSolution* cometPhotometry = nullptr);
	void clearSolution();

	EphemerisTableModel* table() { return &m_table; }

	bool available() const { return m_solution.has_value(); }

	bool busy() const { return m_busy; }

	QString error() const { return m_error; }

	QString summary() const { return m_summary; }

	int epochCount() const { return m_table.count(); }

	int plannedEpochs() const;
	QString formError() const;
	bool startValid() const;
	bool endValid() const;
	bool timeValid() const;

	bool canGenerate() const { return available() && !m_busy && !m_blocked && formError().isEmpty(); }

	void setBlocked(bool value)
	{
		if (m_blocked != value) {
			m_blocked = value;
			emit changed();
		}
	}

	int layout() const { return static_cast<int>(m_table.layout()); }

	void setLayout(int value);

	int center() const { return m_center; }

	void setCenter(int value);

	QString station() const { return m_station; }

	void setStation(const QString& value);

	QString startDate() const { return m_startDate; }

	void setStartDate(const QString& value);

	QString endDate() const { return m_endDate; }

	void setEndDate(const QString& value);

	QString startTime() const { return m_startTime; }

	void setStartTime(const QString& value);

	double stepDays() const { return m_stepDays; }

	void setStepDays(double value);

	int plane() const { return m_plane; }

	void setPlane(int value);

	int scale() const { return m_scale; }

	void setScale(int value);

	bool lightTime() const { return m_lightTime; }

	void setLightTime(bool value);
	QString scaleName() const;

	QStringList exportFormats() const { return {QStringLiteral("TXT"), QStringLiteral("JSON")}; }

	Q_INVOKABLE void generate();
	Q_INVOKABLE void cancel();
	Q_INVOKABLE QString suggestedFileName(int format) const;
	Q_INVOKABLE bool save(const QUrl& file, int format);

   signals:
	void changed();
	void formChanged();

   private:
	void finish(quint64 generation, orbital_core::EphemerisRequest request,
	            std::shared_ptr<orbital_core::EphemerisResult> result);
	void fail(quint64 generation, const QString& message);
	std::optional<orbital_core::EphemerisRequest> request(QString* error) const;
	QString describe() const;
	QByteArray exportText() const;
	QByteArray exportJson() const;

	EphemerisTableModel m_table;
	std::optional<orbital_core::EphemerisSolution> m_solution;
	QString m_designation;
	int m_center{Topocentric};
	QString m_station{QStringLiteral("500")};
	QString m_startDate;
	QString m_endDate;
	QString m_startTime{QStringLiteral("00:00:00")};
	double m_stepDays{1.0};
	int m_plane{0};
	int m_scale{0};
	bool m_lightTime{true};

	bool m_busy{};
	bool m_blocked{};
	QString m_error;
	QString m_summary;
	orbital_core::EphemerisRequest m_lastRequest;
	orbital_core::EphemerisResult m_lastResult;  // flags only; the rows live in m_table
	quint64 m_generation{};
	std::jthread m_worker;
};
