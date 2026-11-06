#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QVariantMap>
#include "orbital_core/Residuals.h"
#include <orbital_core/FitPipeline.h>

// Read-only fit diagnostics and residual charts, owning a snapshot of the
// residuals so QML reads do not race the worker thread.
class DiagnosticsModel : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(double rms READ rms NOTIFY changed)
	Q_PROPERTY(double optimizationTime READ optimizationTime NOTIFY changed)
	Q_PROPERTY(double arcLength READ arcLength NOTIFY changed)
	Q_PROPERTY(QString timeScale READ timeScale NOTIFY changed)

   public:
	explicit DiagnosticsModel(QObject* parent = nullptr);
	~DiagnosticsModel() override;

	void setFromFit(const Residuals& residuals, double optimizationTime, double arcLength);
	void setFromResult(const orbital_core::FitPipelineResult& result, double arcLength);

	void clear();

	[[nodiscard]] double rms() const noexcept { return m_rms; }

	[[nodiscard]] double optimizationTime() const noexcept { return m_optimizationTime; }

	[[nodiscard]] double arcLength() const noexcept { return m_arcLength; }

	QString timeScale() const { return m_timeScale; }

	// Chart fillers (Q_INVOKABLE) — push data into QtCharts QXYSeries via
	// a single bulk replace().  Return aux stats so the QML side doesn't
	// need to re-iterate the cached residuals to compute axis ranges.
	// jdMin/jdMax are MJD in timeScale(); setFromResult uses observation UTC.
	// Chart X values are days from timeOriginMjd (midnight on the first date).
	// Histogram bins have equal widths from zero through histogramMax, with
	// the maximum included in the last bin. Radar/missing rows are omitted.
	Q_INVOKABLE QVariantMap fillResidualCharts(QObject* draSeries, QObject* ddecSeries, QObject* excludedDra = nullptr,
	                                           QObject* excludedDdec = nullptr) const;

   signals:
	void changed();

   private:
	struct ChartSample {
		double mjd;
		double ra;
		double dec;
		bool included;
	};

	void appendSample(const Residual& residual, double jd, bool included);
	std::vector<ChartSample> m_samples;
	QString m_timeScale{QStringLiteral("TDB")};

	double m_rms{0.0};
	double m_optimizationTime{0.0};
	double m_arcLength{0.0};
};
