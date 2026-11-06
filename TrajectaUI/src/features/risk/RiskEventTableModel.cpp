#include "RiskEventTableModel.h"
#include <QDate>
#include <QTime>
#include <orbital_core/TimeScale.h>
#include <cmath>

namespace {
	QString utc(double tdb)
	{
		const double value = orbital_core::time::fromMjdTdb(tdb, orbital_core::time::Scale::UTC);
		if (!std::isfinite(value))
			return QStringLiteral("-");
		auto day = static_cast<qint64>(std::floor(value));
		auto ms = static_cast<qint64>(std::llround((value - day) * 86400000));
		if (ms >= 86400000) {
			++day;
			ms = 0;
		}
		return QDate::fromJulianDay(day + 2400001).toString(Qt::ISODate) + QLatin1Char(' ') +
		       QTime(0, 0).addMSecs(static_cast<int>(ms)).toString(QStringLiteral("HH:mm:ss.zzz"));
	}

	QString number(double value, int digits = 6)
	{
		return std::isfinite(value) ? QString::number(value, 'g', digits) : QStringLiteral("-");
	}
}  // namespace

void RiskEventTableModel::setEvents(const std::vector<orbital_core::hazard::McEvent>& events)
{
	beginResetModel();
	m_rows.clear();
	for (const auto& e : events)
		m_rows.push_back(
		    {{QString::number(m_rows.size() + 1), utc(e.startMjdTdb), utc(e.endMjdTdb), QString::number(e.entryClones),
		      number(e.probability), number(e.probabilityLow), number(e.probabilityHigh)},
		     tr("Physical crossings grouped by time; a gap over 24 hours starts another group. Nominal is not an MC "
		        "trial."),
		     e.entryClones > 0});
	endResetModel();
	emit countChanged();
}

QVariant RiskEventTableModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || index.row() >= count() || index.column() >= columnCount())
		return {};
	const auto& row = m_rows[index.row()];
	if (role == Qt::DisplayRole)
		return row.cells[index.column()];
	if (role == Qt::ToolTipRole)
		return row.detail;
	if (role == Qt::UserRole)
		return row.entry;
	return {};
}

QVariant RiskEventTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole || orientation != Qt::Horizontal || section < 0 || section >= columnCount())
		return {};
	const QStringList labels{tr("Event"), tr("First crossing, UTC"), tr("Last crossing, UTC"), tr("Entries"),
	                         tr("P"),     tr("95% lower"),           tr("95% upper")};
	return labels[section];
}

QHash<int, QByteArray> RiskEventTableModel::roleNames() const
{
	return {{Qt::DisplayRole, "display"}, {Qt::ToolTipRole, "detail"}, {Qt::UserRole, "hasEntry"}};
}
