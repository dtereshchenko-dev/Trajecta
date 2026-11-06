#include "StationResidualModel.h"

#include <QMap>
#include <algorithm>
#include <cmath>
#include <limits>

int StationResidualModel::rowCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : count(); }

QVariant StationResidualModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || index.column() != 0 || index.row() < 0 || index.row() >= count())
		return {};
	const auto& row = m_rows[index.row()];
	switch (role) {
		case StationRole:
			return row.station;
		case IncludedRole:
			return row.included;
		case ExcludedRole:
			return row.excluded;
		case ValidRole:
			return row.valid;
		case RmsRole:
			return row.rms;
		case MeanRaRole:
			return row.meanRa;
		case MeanDecRole:
			return row.meanDec;
		default:
			return {};
	}
}

QHash<int, QByteArray> StationResidualModel::roleNames() const
{
	return {{StationRole, "station"}, {IncludedRole, "included"}, {ExcludedRole, "excluded"}, {ValidRole, "valid"},
	        {RmsRole, "rms"},         {MeanRaRole, "meanRa"},     {MeanDecRole, "meanDec"}};
}

void StationResidualModel::setObservations(std::span<const StationResidualSample> observations)
{
	struct Accumulator {
		Row row;
		double sumRa{};
		double sumDec{};
		double residualNorm{};
	};

	QMap<QString, Accumulator> groups;
	for (const auto& observation : observations) {
		auto& group = groups[observation.station];
		group.row.station = observation.station;
		if (!observation.included) {
			++group.row.excluded;
			continue;
		}
		++group.row.included;
		if (!observation.residual || !std::isfinite(observation.residual->x()) ||
		    !std::isfinite(observation.residual->y()))
			continue;
		++group.row.valid;
		group.sumRa += observation.residual->x();
		group.sumDec += observation.residual->y();
		group.residualNorm = std::hypot(group.residualNorm, observation.residual->x(), observation.residual->y());
	}
	beginResetModel();
	m_rows.clear();
	m_stations = groups.keys();
	m_included = m_excluded = m_missing = 0;
	double maxMean = 0.1;
	for (auto& group : groups) {
		auto& row = group.row;
		m_included += row.included;
		m_excluded += row.excluded;
		m_missing += row.included - row.valid;
		if (row.valid) {
			row.meanRa = group.sumRa / row.valid;
			row.meanDec = group.sumDec / row.valid;
			row.rms = group.residualNorm / std::sqrt(2.0 * row.valid);
			maxMean = std::max({maxMean, std::abs(row.meanRa), std::abs(row.meanDec)});
		} else {
			row.rms = row.meanRa = row.meanDec = std::numeric_limits<double>::quiet_NaN();
		}
		m_rows.push_back(row);
	}
	const double magnitude = std::pow(10.0, std::floor(std::log10(maxMean)));
	for (const double step : {1.0, 2.0, 5.0, 10.0}) {
		m_meanLimit = step * magnitude;
		if (m_meanLimit >= maxMean)
			break;
	}
	sortRows();
	endResetModel();
	emit changed();
}

void StationResidualModel::sortRows()
{
	const auto value = [this](const Row& row) {
		switch (m_sortMode) {
			case Rms:
				return row.valid ? row.rms : -1.0;
			case MeanOffset:
				return row.valid ? std::hypot(row.meanRa, row.meanDec) : -1.0;
			default:
				return static_cast<double>(row.included);
		}
	};
	std::ranges::sort(m_rows, [&](const Row& a, const Row& b) {
		if (m_sortMode != StationCode && value(a) != value(b))
			return value(a) > value(b);
		return a.station < b.station;
	});
}

void StationResidualModel::setSortMode(SortMode mode)
{
	if (mode < IncludedCount || mode > StationCode || mode == m_sortMode)
		return;
	beginResetModel();
	m_sortMode = mode;
	sortRows();
	endResetModel();
	emit sortModeChanged();
}
