#include "EphemerisTableModel.h"

#include <QDate>
#include <QVariantMap>
#include <orbital_core/Constants.h>
#include <cmath>
#include <limits>

namespace {

	constexpr double kRateToArcsecPerHour = 206264.80624709636 / 24.0;
	constexpr double kAuPerDayToKmPerSecond = 149597870.7 / 86400.0;

	struct OptionalColumn {
		const char* key;
		const char* label;
		EphemerisTableModel::Column column;
		bool enabledByDefault;
	};

	constexpr OptionalColumn kOptionalColumns[] = {
	    {"distance", QT_TR_NOOP("Observer distance"), EphemerisTableModel::Distance, true},
	    {"solarDistance", QT_TR_NOOP("Solar distance"), EphemerisTableModel::SolarDistance, true},
	    {"magnitude", QT_TR_NOOP("Magnitude"), EphemerisTableModel::Magnitude, true},
	    {"elongation", QT_TR_NOOP("Solar elongation"), EphemerisTableModel::Elongation, true},
	    {"longitudeRate", QT_TR_NOOP("RA / longitude rate"), EphemerisTableModel::LongitudeRate, false},
	    {"latitudeRate", QT_TR_NOOP("Dec / latitude rate"), EphemerisTableModel::LatitudeRate, false},
	    {"distanceRate", QT_TR_NOOP("Observer range rate"), EphemerisTableModel::DistanceRate, false},
	    {"solarDistanceRate", QT_TR_NOOP("Solar range rate"), EphemerisTableModel::SolarDistanceRate, false},
	    {"lightTime", QT_TR_NOOP("Light time"), EphemerisTableModel::LightTime, false},
	    {"phase", QT_TR_NOOP("Phase angle"), EphemerisTableModel::Phase, false},
	    {"skyMotion", QT_TR_NOOP("Sky motion"), EphemerisTableModel::SkyMotion, false},
	    {"skyPositionAngle", QT_TR_NOOP("Sky motion position angle"), EphemerisTableModel::SkyPositionAngle, false},
	    {"ellipseMajor", QT_TR_NOOP("Error ellipse (3σ)"), EphemerisTableModel::EllipseMajor, false},
	    {"ellipseMinor", QT_TR_NOOP("Error ellipse minor axis (3σ)"), EphemerisTableModel::EllipseMinor, false},
	    {"ellipseAngle", QT_TR_NOOP("Error ellipse orientation"), EphemerisTableModel::EllipseAngle, false}};

	QString scaleName(orbital_core::time::Scale scale)
	{
		switch (scale) {
			case orbital_core::time::Scale::TT:
				return QStringLiteral("TT");
			case orbital_core::time::Scale::TDB:
				return QStringLiteral("TDB");
			case orbital_core::time::Scale::UTC:
			default:
				return QStringLiteral("UTC");
		}
	}

	QString sexagesimal(double value, int decimals, int wholeWidth, bool signedValue, qint64 wrap)
	{
		const qint64 factor = decimals == 3 ? 1000 : decimals == 2 ? 100 : decimals == 1 ? 10 : 1;
		const qint64 perMinute = 60 * factor;
		const qint64 perWhole = 3600 * factor;
		qint64 total = std::llround(std::abs(value) * 3600.0 * static_cast<double>(factor));
		if (wrap > 0)
			total %= wrap * perWhole;
		const qint64 whole = total / perWhole;
		const qint64 minutes = (total % perWhole) / perMinute;
		const qint64 seconds = total % perMinute;
		QString text;
		if (signedValue)
			text += value < 0.0 && total > 0 ? QLatin1Char('-') : QLatin1Char('+');
		text += QStringLiteral("%1 %2 %3")
		            .arg(whole, wholeWidth, 10, QLatin1Char('0'))
		            .arg(minutes, 2, 10, QLatin1Char('0'))
		            .arg(static_cast<double>(seconds) / static_cast<double>(factor), decimals ? 3 + decimals : 2, 'f',
		                 decimals, QLatin1Char('0'));
		return text;
	}

	QString fixed(double value, int decimals)
	{
		return std::isfinite(value) ? QString::number(value, 'f', decimals) : QStringLiteral("—");
	}

	QString significant(double value)
	{
		return std::isfinite(value) ? QString::number(value, 'g', 8) : QStringLiteral("—");
	}

}  // namespace

EphemerisTableModel::EphemerisTableModel(QObject* parent) : QAbstractTableModel(parent)
{
	for (const auto& option : kOptionalColumns)
		m_enabled[0][option.column] = option.enabledByDefault;
}

std::vector<orbital_core::EphemerisRow> EphemerisTableModel::takeRows()
{
	beginResetModel();
	auto rows = std::move(m_rows);
	endResetModel();
	emit changed();
	return rows;
}

int EphemerisTableModel::columns() const
{
	int count = m_layout == Layout::Cartesian ? 7 : 3;
	for (const auto& option : kOptionalColumns)
		count += m_enabled[static_cast<int>(m_layout)][option.column] ? 1 : 0;
	return count;
}

int EphemerisTableModel::logicalColumn(int column) const
{
	if (column < 0 || column >= columns())
		return -1;
	const int baseColumns = m_layout == Layout::Cartesian ? 7 : 3;
	if (m_layout == Layout::Cartesian && column < baseColumns)
		return column == 0 ? Date : X + column - 1;
	if (column < baseColumns)
		return column;
	column -= baseColumns;
	for (const auto& option : kOptionalColumns)
		if (m_enabled[static_cast<int>(m_layout)][option.column] && column-- == 0)
			return option.column;
	return -1;
}

bool EphemerisTableModel::includesColumn(Column column) const
{
	if (column == Date)
		return true;
	if (m_layout == Layout::Cartesian && column >= X && column <= Vz)
		return true;
	if (m_layout == Layout::Spherical && (column == Longitude || column == Latitude))
		return true;
	return column >= Distance && column <= RangeRateUncertainty && m_enabled[static_cast<int>(m_layout)][column];
}

bool EphemerisTableModel::uncertaintySelected() const
{
	for (int column = LongitudeUncertainty; column <= RangeRateUncertainty; ++column)
		if (includesColumn(static_cast<Column>(column)))
			return true;
	return false;
}

QVariantList EphemerisTableModel::columnOptions() const
{
	QVariantList result;
	for (const auto& option : kOptionalColumns) {
		if (option.column == EllipseMinor || option.column == EllipseAngle)
			continue;
		result.append(QVariantMap{{"id", int(option.column)},
		                          {"key", QString::fromLatin1(option.key)},
		                          {"label", tr(option.label)},
		                          {"enabled", m_enabled[static_cast<int>(m_layout)][option.column]}});
	}
	return result;
}

void EphemerisTableModel::setColumnEnabled(int column, bool enabled)
{
	for (const auto& option : kOptionalColumns) {
		if (option.column != column || m_enabled[static_cast<int>(m_layout)][option.column] == enabled)
			continue;
		beginResetModel();
		auto& enabledColumns = m_enabled[static_cast<int>(m_layout)];
		if (option.column >= EllipseMajor && option.column <= EllipseAngle) {
			enabledColumns[EllipseMajor] = enabled;
			enabledColumns[EllipseMinor] = enabled;
			enabledColumns[EllipseAngle] = enabled;
		} else {
			enabledColumns[option.column] = enabled;
		}
		endResetModel();
		emit changed();
		return;
	}
}

int EphemerisTableModel::columnWidth(int column) const
{
	switch (logicalColumn(column)) {
		case Date:
			return 176;
		case Longitude:
			return 132;
		case Latitude:
			return 124;
		case LongitudeRate:
			return 150;
		case Magnitude:
			return 100;
		case Elongation:
		case Phase:
			return 108;
		case SkyPositionAngle:
			return 162;
		case LongitudeUncertainty:
		case LatitudeUncertainty:
		case EllipseMajor:
		case EllipseMinor:
		case EllipseAngle:
		case EllipseArea:
		case SkyUncertainty:
		case RangeUncertainty:
		case RangeRateUncertainty:
			return 165;
		default:
			return 130;
	}
}

QString EphemerisTableModel::calendar(double mjd, int secondDecimals)
{
	const qint64 unitsPerDay = secondDecimals < 0 ? 1440 : secondDecimals == 0 ? 86400 : 86400000;
	double day = std::floor(mjd);
	qint64 units = std::llround((mjd - day) * static_cast<double>(unitsPerDay));
	if (units >= unitsPerDay) {
		units -= unitsPerDay;
		day += 1.0;
	}
	const QDate date = QDate::fromJulianDay(static_cast<qint64>(day) + 2400001);
	QString text = date.toString(QStringLiteral("yyyy-MM-dd"));
	if (secondDecimals < 0)
		return text + QStringLiteral(" %1:%2")
		                  .arg(units / 60, 2, 10, QLatin1Char('0'))
		                  .arg(units % 60, 2, 10, QLatin1Char('0'));
	const qint64 perSecond = unitsPerDay / 86400;
	const qint64 seconds = units / perSecond;
	text += QStringLiteral(" %1:%2:%3")
	            .arg(seconds / 3600, 2, 10, QLatin1Char('0'))
	            .arg((seconds / 60) % 60, 2, 10, QLatin1Char('0'))
	            .arg(seconds % 60, 2, 10, QLatin1Char('0'));
	if (secondDecimals > 0)
		text += QStringLiteral(".%1").arg(units % perSecond, 3, 10, QLatin1Char('0'));
	return text;
}

double EphemerisTableModel::skyMotion(const orbital_core::EphemerisRow& row)
{
	return std::hypot(row.longitudeRate, row.latitudeRate) * kRateToArcsecPerHour;
}

double EphemerisTableModel::skyPositionAngle(const orbital_core::EphemerisRow& row)
{
	if (!std::isfinite(row.longitudeRate) || !std::isfinite(row.latitudeRate) ||
	    (row.longitudeRate == 0.0 && row.latitudeRate == 0.0))
		return std::numeric_limits<double>::quiet_NaN();
	double angle = std::atan2(row.longitudeRate, row.latitudeRate) * 180.0 / Constants::PI;
	return angle < 0.0 ? angle + 360.0 : angle;
}

void EphemerisTableModel::setRows(std::vector<orbital_core::EphemerisRow> rows, orbital_core::EphemerisPlane plane,
                                  orbital_core::time::Scale scale, bool heliocentric,
                                  orbital_core::EphemerisMagnitudeType magnitudeType)
{
	beginResetModel();
	m_rows = std::move(rows);
	m_plane = plane;
	m_scale = scale;
	m_heliocentric = heliocentric;
	m_magnitudeType = magnitudeType;
	m_secondDecimals = -1;
	for (const auto& row : m_rows) {
		const double seconds = (row.mjd - std::floor(row.mjd)) * 86400.0;
		const double rounded = std::round(seconds);
		if (std::abs(seconds - rounded) > 1.0e-3) {
			m_secondDecimals = 3;
			break;
		}
		if (std::fmod(rounded, 60.0) != 0.0)
			m_secondDecimals = 0;
	}
	endResetModel();
	emit headerDataChanged(Qt::Horizontal, 0, columns() - 1);
	emit changed();
}

void EphemerisTableModel::setLayout(Layout layout)
{
	if (m_layout == layout)
		return;
	beginResetModel();
	m_layout = layout;
	endResetModel();
	emit headerDataChanged(Qt::Horizontal, 0, columns() - 1);
	emit changed();
}

void EphemerisTableModel::clear()
{
	beginResetModel();
	m_rows.clear();
	endResetModel();
	emit changed();
}

QString EphemerisTableModel::header(int column) const
{
	const int logical = logicalColumn(column);
	return logical >= Date && logical < ColumnCount ? headerForColumn(static_cast<Column>(logical)) : QString();
}

QString EphemerisTableModel::headerForColumn(Column logical) const
{
	if (logical == Date)
		return tr("Date %1").arg(scaleName(m_scale));
	if (logical >= X && logical <= Vz) {
		static const char* names[] = {"X, AU", "Y, AU", "Z, AU", "Vx, AU/d", "Vy, AU/d", "Vz, AU/d"};
		return logical >= X && logical <= Vz ? tr(names[logical - X]) : QString();
	}
	const bool ecliptic = m_plane == orbital_core::EphemerisPlane::Ecliptic;
	switch (logical) {
		case Longitude:
			return ecliptic ? tr("λ, ° ′ ″") : tr("RA, h m s");
		case Latitude:
			return ecliptic ? tr("β, ° ′ ″") : tr("Dec, ° ′ ″");
		case Distance:
			return tr("Δ, AU");
		case SolarDistance:
			return tr("r, AU");
		case Magnitude:
			if (m_magnitudeType == orbital_core::EphemerisMagnitudeType::CometTotal)
				return tr("Total, mag");
			if (m_magnitudeType == orbital_core::EphemerisMagnitudeType::CometNuclear)
				return tr("Nuclear, mag");
			return tr("V, mag");
		case Elongation:
			return tr("Elongation, °");
		case LongitudeRate:
			return ecliptic ? tr("dλ cos β, ″/h") : tr("dRA cos Dec, ″/h");
		case LatitudeRate:
			return ecliptic ? tr("dβ, ″/h") : tr("dDec, ″/h");
		case DistanceRate:
			return tr("Δ rate, km/s");
		case SolarDistanceRate:
			return tr("r rate, km/s");
		case LightTime:
			return tr("Light time, min");
		case Phase:
			return tr("Phase, °");
		case SkyMotion:
			return tr("Sky motion, ″/h");
		case SkyPositionAngle:
			return tr("Sky motion PA, °");
		case LongitudeUncertainty:
			return ecliptic ? tr("λ cos β 3σ, ″") : tr("RA cos Dec 3σ, ″");
		case LatitudeUncertainty:
			return ecliptic ? tr("β 3σ, ″") : tr("Dec 3σ, ″");
		case EllipseMajor:
			return tr("Ellipse major 3σ, ″");
		case EllipseMinor:
			return tr("Ellipse minor 3σ, ″");
		case EllipseAngle:
			return tr("Ellipse angle, °");
		case EllipseArea:
			return tr("Ellipse area 3σ, ″²");
		case SkyUncertainty:
			return tr("Sky 3σ RSS, ″");
		case RangeUncertainty:
			return tr("Range 3σ, km");
		case RangeRateUncertainty:
			return tr("Range rate 3σ, km/s");
		default:
			return {};
	}
}

QString EphemerisTableModel::cell(int row, int column) const
{
	const int logical = logicalColumn(column);
	return logical >= Date && logical < ColumnCount ? cellForColumn(row, static_cast<Column>(logical)) : QString();
}

QString EphemerisTableModel::cellForColumn(int row, Column logical) const
{
	if (row < 0 || row >= count() || logical < Date || logical >= ColumnCount)
		return {};
	const auto& e = m_rows[static_cast<std::size_t>(row)];
	if (logical == Date)
		return calendar(e.mjd, m_secondDecimals);
	if (logical >= X && logical <= Vz) {
		const double values[] = {e.position.x, e.position.y, e.position.z, e.velocity.x, e.velocity.y, e.velocity.z};
		return fixed(values[logical - X], 10);
	}
	const bool ecliptic = m_plane == orbital_core::EphemerisPlane::Ecliptic;
	switch (logical) {
		case Longitude:
			return ecliptic ? sexagesimal(e.longitude * 180.0 / Constants::PI, 2, 3, false, 360)
			                : sexagesimal(e.longitude * 12.0 / Constants::PI, 3, 2, false, 24);
		case Latitude:
			return sexagesimal(e.latitude * 180.0 / Constants::PI, 2, 2, true, 0);
		case Distance:
			return fixed(e.distanceAu, 8);
		case SolarDistance:
			return fixed(e.heliocentricAu, 6);
		case Magnitude:
			return e.visualMagnitude ? fixed(*e.visualMagnitude, 2) : QStringLiteral("—");
		case Elongation:
			return m_heliocentric ? QStringLiteral("—") : fixed(e.elongation * 180.0 / Constants::PI, 3);
		case LongitudeRate:
			return fixed(e.longitudeRate * kRateToArcsecPerHour, 3);
		case LatitudeRate:
			return fixed(e.latitudeRate * kRateToArcsecPerHour, 3);
		case DistanceRate:
			return fixed(e.distanceRateAuPerDay * kAuPerDayToKmPerSecond, 6);
		case SolarDistanceRate:
			return fixed(e.heliocentricRateAuPerDay * kAuPerDayToKmPerSecond, 6);
		case LightTime:
			return fixed(e.oneWayLightTimeDays * 1440.0, 6);
		case Phase:
			return m_heliocentric ? QStringLiteral("—") : fixed(e.phase * 180.0 / Constants::PI, 3);
		case SkyMotion:
			return fixed(skyMotion(e), 3);
		case SkyPositionAngle:
			return fixed(skyPositionAngle(e), 3);
		case LongitudeUncertainty:
			return e.uncertainty ? significant(e.uncertainty->longitude3SigmaArcsec) : QStringLiteral("—");
		case LatitudeUncertainty:
			return e.uncertainty ? significant(e.uncertainty->latitude3SigmaArcsec) : QStringLiteral("—");
		case EllipseMajor:
			return e.uncertainty ? significant(e.uncertainty->major3SigmaArcsec) : QStringLiteral("—");
		case EllipseMinor:
			return e.uncertainty ? significant(e.uncertainty->minor3SigmaArcsec) : QStringLiteral("—");
		case EllipseAngle:
			return e.uncertainty ? fixed(e.uncertainty->positionAngleDeg, 3) : QStringLiteral("—");
		case EllipseArea:
			return e.uncertainty ? significant(e.uncertainty->area3SigmaArcsec2) : QStringLiteral("—");
		case SkyUncertainty:
			return e.uncertainty ? significant(e.uncertainty->sky3SigmaArcsec) : QStringLiteral("—");
		case RangeUncertainty:
			return e.uncertainty ? significant(e.uncertainty->range3SigmaKm) : QStringLiteral("—");
		case RangeRateUncertainty:
			return e.uncertainty ? significant(e.uncertainty->rangeRate3SigmaKmPerSecond) : QStringLiteral("—");
		default:
			return {};
	}
}

int EphemerisTableModel::rowCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : count(); }

int EphemerisTableModel::columnCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : columns(); }

QVariant EphemerisTableModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || role != Qt::DisplayRole)
		return {};
	return cell(index.row(), index.column());
}

QVariant EphemerisTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (orientation != Qt::Horizontal || role != Qt::DisplayRole)
		return {};
	return header(section);
}
