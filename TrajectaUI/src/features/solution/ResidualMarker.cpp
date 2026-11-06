#include "ResidualMarker.h"

#include <QtCharts/QChart>
#include <QtCharts/QScatterSeries>
#include <cmath>

void ResidualMarker::setSeries(QXYSeries* series)
{
	if (m_series == series)
		return;
	if (m_series)
		m_series->setPointsVisible(false);
	m_series = series;
	if (m_series) {
		m_series->setMarkerSize(3.2);
		m_series->setPointsVisible(true);
	}
	emit seriesChanged();
}

int ResidualMarker::pointIndexAt(QXYSeries* series, const QPointF& chartPosition) const
{
	if (!series || !series->isVisible() || !series->chart() || !series->chart()->plotArea().contains(chartPosition))
		return -1;
	const auto value = series->chart()->mapToValue(chartPosition, series);
	const auto unitsPerPixel = series->chart()->mapToValue(chartPosition + QPointF(1, 1), series) - value;
	if (!std::isfinite(unitsPerPixel.x()) || !std::isfinite(unitsPerPixel.y()) || unitsPerPixel.x() == 0 ||
	    unitsPerPixel.y() == 0)
		return -1;
	const double radius = series->markerSize() / (qobject_cast<QScatterSeries*>(series) ? 2 : 1);
	double nearestSquared = (radius + 1) * (radius + 1);
	int nearest = -1;
	const auto points = series->points();
	for (int i = 0; i < points.size(); ++i) {
		const double dx = (points[i].x() - value.x()) / unitsPerPixel.x();
		const double dy = (points[i].y() - value.y()) / unitsPerPixel.y();
		const double distanceSquared = dx * dx + dy * dy;
		if (distanceSquared <= nearestSquared) {
			nearestSquared = distanceSquared;
			nearest = i;
		}
	}
	return nearest;
}
