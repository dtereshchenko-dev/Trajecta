#pragma once

#include <QPointer>
#include <QQmlEngine>
#include <QtCharts/QXYSeries>

// Hit testing shares the chart's coordinate transform and ignores invisible
// connecting lines. Qt 6.9 LineSeries vector points use markerSize as a radius;
// ScatterSeries selection overlays use markerSize as a diameter.
class ResidualMarker : public QObject {
	Q_OBJECT
	QML_ELEMENT
	Q_PROPERTY(QXYSeries* series READ series WRITE setSeries NOTIFY seriesChanged REQUIRED)

   public:
	using QObject::QObject;

	QXYSeries* series() const { return m_series; }

	void setSeries(QXYSeries* series);
	Q_INVOKABLE int pointIndexAt(QXYSeries* series, const QPointF& chartPosition) const;

   signals:
	void seriesChanged();

   private:
	QPointer<QXYSeries> m_series;
};
