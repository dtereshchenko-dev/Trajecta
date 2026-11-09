#include "CloseApproachTableModel.h"
#include <orbital_core/Constants.h>
#include <QDate>
#include <QLocale>
#include <QTime>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace {
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	const QStringList distances{"au", "km", "LD", "mi", "Re"};
	const QStringList velocities{"km/s", "mi/h", "au/d"};

	QString number(double value)
	{
		return std::isfinite(value) ? QString::number(value, 'g', 9) : QStringLiteral("\u2014");
	}

	QString epoch(double mjd)
	{
		if (!std::isfinite(mjd))
			return QStringLiteral("\u2014");
		const auto totalMinutes = static_cast<qint64>(std::llround(mjd * 1440.0));
		const auto day = static_cast<qint64>(std::floor(double(totalMinutes) / 1440.0));
		const auto date = QDate::fromJulianDay(day + 2400001);
		return QLocale::c().toString(date, QStringLiteral("yyyy-MMM-dd")) + QLatin1Char(' ') +
		       QTime(0, 0).addSecs(static_cast<int>((totalMinutes - day * 1440) * 60)).toString("HH:mm");
	}

	QString uncertainty(double days)
	{
		if (!std::isfinite(days) || days < 0.0)
			return QStringLiteral("\u2014");
		if (days * 1440.0 < 1.0)
			return QStringLiteral("<00:01");
		if (days > 1000000000.0)
			return QString::number(days, 'g', 5) + QStringLiteral(" d");
		const auto minutes = static_cast<qint64>(std::ceil(days * 1440.0));
		return (minutes >= 1440 ? QString::number(minutes / 1440) + QLatin1Char('_') : QString()) +
		       QStringLiteral("%1:%2")
		           .arg((minutes / 60) % 24, 2, 10, QLatin1Char('0'))
		           .arg(minutes % 60, 2, 10, QLatin1Char('0'));
	}

	double numericValue(const orbital_core::hazard::NominalCloseApproach& row, int col)
	{
		const double nan = std::numeric_limits<double>::quiet_NaN();
		switch (col) {
			case 0:
				return row.mjdTdb;
			case 2:
				return row.encounter.distAu;
			case 3:
				return row.minimumDistanceAu.value_or(nan);
			case 4:
				return row.maximumDistanceAu.value_or(nan);
			case 5:
				return row.relativeSpeedKmS;
			case 6:
				return row.encounter.vInfKmS > 0.0 ? row.encounter.vInfKmS : nan;
			default:
				return nan;
		}
	}
}  // namespace

QHash<int, QByteArray> CloseApproachTableModel::roleNames() const
{
	return {{Qt::DisplayRole, "display"}};
}

QVariant CloseApproachTableModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || index.row() >= count() || index.column() >= 7)
		return {};
	if (role == Qt::DisplayRole)
		return cell(index.row(), index.column());
	return {};
}

QVariant CloseApproachTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	return orientation == Qt::Horizontal && role == Qt::DisplayRole && section >= 0 && section < 7
	           ? QVariant(header(section))
	           : QVariant();
}

QString CloseApproachTableModel::header(int column) const
{
	const QStringList headers{tr("Date/Time (TDB)"),
	                          tr("Body"),
	                          tr("Nominal distance\n(%1)").arg(distances[m_distanceUnit]),
	                          tr("Minimum distance\n(%1)").arg(distances[m_distanceUnit]),
	                          tr("Maximum distance\n(%1)").arg(distances[m_distanceUnit]),
	                          tr("V-relative\n(%1)").arg(velocities[m_velocityUnit]),
	                          tr("V-infinity\n(%1)").arg(velocities[m_velocityUnit])};
	return headers.value(column);
}

QString CloseApproachTableModel::cell(int index, int column) const
{
	if (index < 0 || index >= count())
		return {};
	const auto& row = m_rows[index];
	if (column == 0)
		return epoch(row.mjdTdb) + (row.timeUncertaintyDays
		                                ? QStringLiteral(" \u00b1 ") + uncertainty(*row.timeUncertaintyDays)
		                                : QString());
	if (column == 1)
		return QString::fromStdString(std::string(orbital_core::hazard::longLabel(row.encounter.planet)));
	const double distanceFactors[]{1.0, auKm, auKm / Constants::LUNAR_DISTANCE_KM, auKm / 1.609344,
	                               auKm / Constants::EARTH_EQ_RADIUS_KM};
	const double speedFactors[]{1.0, 3600.0 / 1.609344, 86400.0 / auKm};
	return number(numericValue(row, column) *
	              (column < 5 ? distanceFactors[m_distanceUnit] : speedFactors[m_velocityUnit]));
}

QString CloseApproachTableModel::toPlainText() const
{
	std::array<qsizetype, 7> widths{};
	for (int c = 0; c < columnCount(); ++c) {
		widths[c] = header(c).size();
		for (int row = 0; row < count(); ++row)
			widths[c] = std::max(widths[c], cell(row, c).size());
	}
	QString text;
	for (int row = -1; row < count(); ++row) {
		for (int c = 0; c < columnCount(); ++c) {
			if (c)
				text += "  ";
			const auto value = row < 0 ? header(c).replace('\n', ' ') : cell(row, c);
			text += c < 2 ? value.leftJustified(widths[c]) : value.rightJustified(widths[c]);
		}
		text += '\n';
	}
	return text;
}

void CloseApproachTableModel::sortRows()
{
	std::stable_sort(m_rows.begin(), m_rows.end(), [this](const auto& a, const auto& b) {
		if (m_sortColumn == 1) {
			const auto aa = orbital_core::hazard::longLabel(a.encounter.planet);
			const auto bb = orbital_core::hazard::longLabel(b.encounter.planet);
			return m_descending ? aa > bb : aa < bb;
		}
		const double aa = numericValue(a, m_sortColumn), bb = numericValue(b, m_sortColumn);
		if (!std::isfinite(aa))
			return false;
		if (!std::isfinite(bb))
			return true;
		return m_descending ? aa > bb : aa < bb;
	});
}

void CloseApproachTableModel::setRows(std::vector<orbital_core::hazard::NominalCloseApproach> rows)
{
	beginResetModel();
	m_rows = std::move(rows);
	sortRows();
	endResetModel();
	emit changed();
}

void CloseApproachTableModel::toggleSort(int column)
{
	if (column < 0 || column >= 7)
		return;
	beginResetModel();
	m_descending = column == m_sortColumn ? !m_descending : false;
	m_sortColumn = column;
	sortRows();
	endResetModel();
	emit changed();
}

void CloseApproachTableModel::reformat()
{
	emit headerDataChanged(Qt::Horizontal, 0, 6);
	if (count())
		emit dataChanged(index(0, 0), index(count() - 1, 6));
	emit changed();
}

void CloseApproachTableModel::setDistanceUnit(int value)
{
	if (value < 0 || value >= distances.size() || value == m_distanceUnit)
		return;
	m_distanceUnit = value;
	reformat();
}

void CloseApproachTableModel::setVelocityUnit(int value)
{
	if (value < 0 || value >= velocities.size() || value == m_velocityUnit)
		return;
	m_velocityUnit = value;
	reformat();
}
