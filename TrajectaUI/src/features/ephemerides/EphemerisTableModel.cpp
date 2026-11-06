#include "EphemerisTableModel.h"

#include <QDate>
#include <cmath>

namespace {

	constexpr double kPi = 3.14159265358979323846;

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

}  // namespace

EphemerisTableModel::EphemerisTableModel(QObject* parent) : QAbstractTableModel(parent) {}

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
	if (column == 0)
		return tr("Date %1").arg(scaleName(m_scale));
	if (m_layout == Layout::Cartesian) {
		static const char* names[] = {"", "X, AU", "Y, AU", "Z, AU", "Vx, AU/d", "Vy, AU/d", "Vz, AU/d"};
		return column > 0 && column < 7 ? tr(names[column]) : QString();
	}
	const bool ecliptic = m_plane == orbital_core::EphemerisPlane::Ecliptic;
	switch (column) {
		case 1:
			return ecliptic ? tr("λ, ° ′ ″") : tr("RA, h m s");
		case 2:
			return ecliptic ? tr("β, ° ′ ″") : tr("Dec, ° ′ ″");
		case 3:
			return tr("Δ, AU");
		case 4:
			return tr("r, AU");
		case 5:
			if (m_magnitudeType == orbital_core::EphemerisMagnitudeType::CometTotal)
				return tr("Total, mag");
			if (m_magnitudeType == orbital_core::EphemerisMagnitudeType::CometNuclear)
				return tr("Nuclear, mag");
			return tr("V, mag");
		case 6:
			return tr("Elongation, °");
		default:
			return {};
	}
}

QString EphemerisTableModel::cell(int row, int column) const
{
	if (row < 0 || row >= count() || column < 0 || column >= columns())
		return {};
	const auto& e = m_rows[static_cast<std::size_t>(row)];
	if (column == 0)
		return calendar(e.mjd, m_secondDecimals);
	if (m_layout == Layout::Cartesian) {
		const double values[] = {0.0,          e.position.x, e.position.y, e.position.z,
		                         e.velocity.x, e.velocity.y, e.velocity.z};
		return fixed(values[column], 10);
	}
	const bool ecliptic = m_plane == orbital_core::EphemerisPlane::Ecliptic;
	switch (column) {
		case 1:
			return ecliptic ? sexagesimal(e.longitude * 180.0 / kPi, 2, 3, false, 360)
			                : sexagesimal(e.longitude * 12.0 / kPi, 3, 2, false, 24);
		case 2:
			return sexagesimal(e.latitude * 180.0 / kPi, 2, 2, true, 0);
		case 3:
			return fixed(e.distanceAu, 8);
		case 4:
			return fixed(e.heliocentricAu, 6);
		case 5:
			return e.visualMagnitude ? fixed(*e.visualMagnitude, 2) : QStringLiteral("—");
		case 6:
			return m_heliocentric ? QStringLiteral("—") : fixed(e.elongation * 180.0 / kPi, 3);
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
