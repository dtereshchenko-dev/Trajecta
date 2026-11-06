#include "ObservationReviewModel.h"

#include "util/FormatUtils.h"

#include <QDate>
#include <QDateTime>
#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace {
	enum Roles {
		GroupRole = Qt::UserRole + 1,
		StationRole,
		TimeRole,
		CoordinatesRole,
		ReasonRole,
		IncludedRole,
		CanChangeRole
	};

	QString epoch(const ObsNode& observation)
	{
		QDateTime date;
		if (observation.hasSplitUtc) {
			const auto day = QDate::fromJulianDay(static_cast<qint64>(observation.utcMjdDay) + 2400001);
			date = QDateTime(day, QTime(0, 0), Qt::UTC).addMSecs(qRound64(observation.utcSecondsOfDay * 1000.0));
		} else {
			date = QDateTime::fromMSecsSinceEpoch(qRound64((observation.JD - 2440587.5) * 86400000.0), Qt::UTC);
		}
		const char* scale = observation.timeScale == orbital_core::time::Scale::UT1   ? "UT1"
		                    : observation.timeScale == orbital_core::time::Scale::TT  ? "TT"
		                    : observation.timeScale == orbital_core::time::Scale::TDB ? "TDB"
		                                                                              : "UTC";
		return date.toString("yyyy-MM-dd HH:mm:ss.zzz") + ' ' + QString::fromLatin1(scale);
	}
}  // namespace

ObservationReviewModel::ObservationReviewModel(QObject* parent) : QAbstractListModel(parent) {}

void ObservationReviewModel::setObservations(const orbital_core::ObservationSet& observations)
{
	beginResetModel();
	m_rows.clear();
	m_filter = 0;
	m_initialOpticalCount = static_cast<int>(std::ranges::count_if(
	    observations.all, [](const ObsNode& o) { return !o.radar && o.note2 != 'X' && o.note2 != 'x'; }));
	std::unordered_map<orbital_core::ObservationId, const ObsNode*> byId;
	for (const auto& o : observations.all)
		byId.emplace(o.id, &o);
	std::unordered_map<orbital_core::ObservationId, std::size_t> rowById;
	const auto ensureRow = [&](const orbital_core::ObservationId& id) -> Row* {
		if (const auto existing = rowById.find(id); existing != rowById.end())
			return &m_rows[existing->second];
		const auto found = byId.find(id);
		if (found == byId.end())
			return nullptr;
		const auto& o = *found->second;
		const bool replaced = o.note2 == 'X' || o.note2 == 'x';
		rowById.emplace(id, m_rows.size());
		m_rows.push_back({id,
		                  QString::fromStdString(o.obsCode),
		                  epoch(o),
		                  FormatUtils::ra(o.submittedRaRad()) + "  " + FormatUtils::dec(o.submittedDecRad()),
		                  {},
		                  {},
		                  0,
		                  !replaced,
		                  replaced});
		return &m_rows.back();
	};
	const auto addReason = [](Row& row, orbital_core::ObservationIssue issue, const QString& reason) {
		if (std::ranges::find(row.issues, issue) != row.issues.end())
			return;
		row.issues.push_back(issue);
		if (!row.reason.isEmpty())
			row.reason += "; ";
		row.reason += reason;
	};
	int groupNumber = 0;
	for (const auto& group : observations.coincidentObservations) {
		++groupNumber;
		for (const auto& id : group.observations) {
			auto* row = ensureRow(id);
			if (!row)
				continue;
			row->group = groupNumber;
			for (const auto issue : group.issues) {
				if (issue == orbital_core::ObservationIssue::CoincidentPosition)
					addReason(*row, issue, tr("Same time and position"));
				else if (issue == orbital_core::ObservationIssue::NearbyObservationConflict)
					addReason(*row, issue,
					          tr("Within 3 seconds; conflicting positions; maximum separation from first: %1 arcsec")
					              .arg(group.maxSeparationArcsec, 0, 'g', 5));
			}
		}
	}
	for (const auto& diagnostic : observations.visibilityDiagnostics) {
		auto* row = ensureRow(diagnostic.observation);
		if (!row)
			continue;
		for (const auto issue : diagnostic.issues) {
			switch (issue) {
				case orbital_core::ObservationIssue::BelowHorizon:
					addReason(*row, issue,
					          tr("Object below horizon at submitted time: %1 deg")
					              .arg(*diagnostic.objectAltitudeDeg, 0, 'f', 1));
					break;
				case orbital_core::ObservationIssue::Daylight:
					addReason(
					    *row, issue,
					    tr("Sun above horizon at submitted time: %1 deg").arg(*diagnostic.sunAltitudeDeg, 0, 'f', 1));
					break;
				default:
					break;
			}
		}
	}
	rebuildFilter();
	endResetModel();
	emit filterChanged();
	emit changed();
}

void ObservationReviewModel::clear()
{
	beginResetModel();
	m_rows.clear();
	m_visible.clear();
	m_filter = 0;
	m_initialOpticalCount = 0;
	endResetModel();
	emit filterChanged();
	emit changed();
}

int ObservationReviewModel::rowCount(const QModelIndex& parent) const
{
	return parent.isValid() ? 0 : static_cast<int>(m_visible.size());
}

QVariant ObservationReviewModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || index.row() < 0 || index.row() >= static_cast<int>(m_visible.size()))
		return {};
	const auto& row = m_rows[m_visible[index.row()]];
	switch (role) {
		case GroupRole:
			return row.group;
		case StationRole:
			return row.station;
		case TimeRole:
			return row.time;
		case CoordinatesRole:
			return row.coordinates;
		case ReasonRole:
			return row.reason;
		case IncludedRole:
			return row.included;
		case CanChangeRole:
			return !row.replaced;
		default:
			return {};
	}
}

QHash<int, QByteArray> ObservationReviewModel::roleNames() const
{
	return {{GroupRole, "groupNumber"},       {StationRole, "station"}, {TimeRole, "observationTime"},
	        {CoordinatesRole, "coordinates"}, {ReasonRole, "reason"},   {IncludedRole, "included"},
	        {CanChangeRole, "canChange"}};
}

void ObservationReviewModel::rebuildFilter()
{
	m_visible.clear();
	for (std::size_t i = 0; i < m_rows.size(); ++i) {
		const auto issue = static_cast<orbital_core::ObservationIssue>(m_filter - 1);
		if (m_filter == 0 || std::ranges::find(m_rows[i].issues, issue) != m_rows[i].issues.end())
			m_visible.push_back(i);
	}
}

void ObservationReviewModel::setFilter(int value)
{
	if (value != 0 && value != positionFilter() && value != differentPositionFilter() &&
	    value != belowHorizonFilter() && value != daylightFilter())
		return;
	if (m_filter == value)
		return;
	beginResetModel();
	m_filter = value;
	rebuildFilter();
	endResetModel();
	emit filterChanged();
}

int ObservationReviewModel::excludedCount() const
{
	return static_cast<int>(
	    std::ranges::count_if(m_rows, [](const Row& row) { return !row.included && !row.replaced; }));
}

void ObservationReviewModel::setIncluded(int visibleRow, bool included)
{
	if (visibleRow < 0 || visibleRow >= static_cast<int>(m_visible.size()))
		return;
	auto& row = m_rows[m_visible[visibleRow]];
	if (row.replaced || row.included == included)
		return;
	row.included = included;
	emit dataChanged(index(visibleRow), index(visibleRow), {IncludedRole});
	emit changed();
}

void ObservationReviewModel::setAllShownIncluded(bool included)
{
	bool changedAny = false;
	for (const auto index : m_visible) {
		auto& row = m_rows[index];
		if (!row.replaced && row.included != included) {
			row.included = included;
			changedAny = true;
		}
	}
	if (changedAny) {
		emit dataChanged(index(0), index(static_cast<int>(m_visible.size()) - 1), {IncludedRole});
		emit changed();
	}
}

std::vector<orbital_core::ObservationOverride> ObservationReviewModel::overrides() const
{
	std::vector<orbital_core::ObservationOverride> result;
	for (const auto& row : m_rows)
		if (!row.included && !row.replaced)
			result.push_back({row.id, orbital_core::ObservationDisposition::Exclude});
	return result;
}
