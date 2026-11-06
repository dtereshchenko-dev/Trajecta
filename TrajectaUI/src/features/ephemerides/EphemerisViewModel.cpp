#include "EphemerisViewModel.h"
#include "util/CalculationTime.h"

#include <QDate>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTime>
#include <orbital_core/Types.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace {

	constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;
	constexpr double kRadPerDayToArcsecPerHour = 206264.80624709636 / 24.0;

	orbital_core::time::Scale toScale(int value)
	{
		return value == 1   ? orbital_core::time::Scale::TT
		       : value == 2 ? orbital_core::time::Scale::TDB
		                    : orbital_core::time::Scale::UTC;
	}

	QString scaleLabel(orbital_core::time::Scale scale)
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

	QTime parseTime(const QString& time)
	{
		for (const auto* format : {"HH:mm:ss.zzz", "HH:mm:ss", "HH:mm"}) {
			const QTime clock = QTime::fromString(time.trimmed(), QLatin1String(format));
			if (clock.isValid())
				return clock;
		}
		return {};
	}

	std::optional<double> parseMjd(const QString& date, const QString& time)
	{
		const QDate day = QDate::fromString(date.trimmed(), Qt::ISODate);
		const QTime clock = parseTime(time);
		if (!day.isValid() || !clock.isValid())
			return std::nullopt;
		return static_cast<double>(day.toJulianDay() - 2400001) + clock.msecsSinceStartOfDay() / 86400000.0;
	}

	QString isoDate(double mjd)
	{
		return QDate::fromJulianDay(static_cast<qint64>(std::floor(mjd)) + 2400001).toString(Qt::ISODate);
	}

	QString number(double value) { return std::isfinite(value) ? QString::number(value, 'g', 17) : QString(); }

	QJsonValue jsonNumber(double value)
	{
		return std::isfinite(value) ? QJsonValue(value) : QJsonValue(QJsonValue::Null);
	}

	struct NumericColumn {
		const char* key;
		double (*value)(const orbital_core::EphemerisRow&);
	};

	const std::vector<NumericColumn>& numericColumns()
	{
		using Row = orbital_core::EphemerisRow;
		static const std::vector<NumericColumn> columns{
		    {"mjd_utc", [](const Row& r) { return r.mjdUtc; }},
		    {"mjd_tt", [](const Row& r) { return r.mjdTt; }},
		    {"mjd_tdb", [](const Row& r) { return r.mjdTdb; }},
		    {"lon_deg", [](const Row& r) { return r.longitude * kRadToDeg; }},
		    {"lat_deg", [](const Row& r) { return r.latitude * kRadToDeg; }},
		    {"lon_rate_coslat_arcsec_h", [](const Row& r) { return r.longitudeRate * kRadPerDayToArcsecPerHour; }},
		    {"lat_rate_arcsec_h", [](const Row& r) { return r.latitudeRate * kRadPerDayToArcsecPerHour; }},
		    {"delta_au", [](const Row& r) { return r.distanceAu; }},
		    {"r_au", [](const Row& r) { return r.heliocentricAu; }},
		    {"elongation_deg", [](const Row& r) { return r.elongation * kRadToDeg; }},
		    {"phase_deg", [](const Row& r) { return r.phase * kRadToDeg; }},
		    {"v_mag",
		     [](const Row& r) { return r.visualMagnitude.value_or(std::numeric_limits<double>::quiet_NaN()); }},
		    {"x_au", [](const Row& r) { return r.position.x; }},
		    {"y_au", [](const Row& r) { return r.position.y; }},
		    {"z_au", [](const Row& r) { return r.position.z; }},
		    {"vx_au_d", [](const Row& r) { return r.velocity.x; }},
		    {"vy_au_d", [](const Row& r) { return r.velocity.y; }},
		    {"vz_au_d", [](const Row& r) { return r.velocity.z; }},
		};
		return columns;
	}

	QString columnKey(const NumericColumn& column, orbital_core::EphemerisPlane plane)
	{
		QString key = QString::fromLatin1(column.key);
		if (plane == orbital_core::EphemerisPlane::Equatorial)
			key.replace(QStringLiteral("lon_rate_coslat"), QStringLiteral("ra_rate_cosdec"))
			    .replace(QStringLiteral("lat_rate"), QStringLiteral("dec_rate"))
			    .replace(QStringLiteral("lon_deg"), QStringLiteral("ra_deg"))
			    .replace(QStringLiteral("lat_deg"), QStringLiteral("dec_deg"));
		return key;
	}

}  // namespace

EphemerisViewModel::EphemerisViewModel(QObject* parent) : QObject(parent)
{
	connect(this, &EphemerisViewModel::formChanged, this, &EphemerisViewModel::changed);
	connect(&m_table, &EphemerisTableModel::changed, this, &EphemerisViewModel::changed);
}

EphemerisViewModel::~EphemerisViewModel()
{
	m_worker.request_stop();
	if (m_worker.joinable())
		m_worker.join();
}

void EphemerisViewModel::setSolution(const orbital_core::FitPipelineResult& result, const QString& designation,
                                     std::optional<double> absoluteMagnitude,
                                     const orbital_core::photometry::CometSolution* cometPhotometry)
{
	const auto& orbit = result.fit.fitEpochOrbit ? result.fit.fitEpochOrbit : result.fit.orbit;
	if (result.fit.cancelled || !orbit) {
		clearSolution();
		return;
	}
	if (m_busy) {
		m_worker.request_stop();
		++m_generation;
	}
	orbital_core::EphemerisSolution solution;
	solution.position = orbit->getR();
	solution.velocity = orbit->getV() * Constants::GAUSS_K;
	solution.epochMjdTdb = orbit->epochMjd();
	solution.forces = result.effectiveConfig.nbody;
	solution.nonGrav = result.fit.finalNonGrav;
	solution.slopeG = result.effectiveConfig.photometry.asteroidG;
	if (absoluteMagnitude && std::isfinite(*absoluteMagnitude))
		solution.absoluteMagnitude = absoluteMagnitude;
	if (cometPhotometry) {
		const auto& component =
		    cometPhotometry->total.absoluteMagnitude ? cometPhotometry->total : cometPhotometry->nuclear;
		solution.magnitudeType = cometPhotometry->total.absoluteMagnitude
		                             ? orbital_core::EphemerisMagnitudeType::CometTotal
		                             : orbital_core::EphemerisMagnitudeType::CometNuclear;
		solution.absoluteMagnitude = component.absoluteMagnitude;
		solution.cometSlope = component.slope;
		solution.cometPhaseCoefficient = component.phaseCoefficient;
	}
	m_solution = solution;

	if (designation != m_designation || m_startDate.isEmpty()) {
		const auto startUtc = CalculationTime::defaultStartUtc(result.fit);
		m_startDate = startUtc ? isoDate(*startUtc) : QString();
		m_endDate = startUtc ? isoDate(std::floor(*startUtc) + 14.0) : QString();
	}
	m_designation = designation;
	m_error.clear();
	m_summary.clear();
	m_table.clear();
	emit formChanged();
}

void EphemerisViewModel::clearSolution()
{
	if (m_busy) {
		m_worker.request_stop();
		++m_generation;
	}
	m_solution.reset();
	m_error.clear();
	m_summary.clear();
	m_table.clear();
	emit changed();
}

void EphemerisViewModel::setLayout(int value)
{
	const auto next = value == 1 ? EphemerisTableModel::Layout::Cartesian : EphemerisTableModel::Layout::Spherical;
	if (next == m_table.layout())
		return;
	m_table.setLayout(next);
	emit formChanged();
}

void EphemerisViewModel::setCenter(int value)
{
	value = std::clamp(value, 0, 2);
	if (m_center == value)
		return;
	m_center = value;
	emit formChanged();
}

void EphemerisViewModel::setStation(const QString& value)
{
	const QString code = value.trimmed().toUpper();
	if (m_station == code)
		return;
	m_station = code;
	emit formChanged();
}

void EphemerisViewModel::setStartDate(const QString& value)
{
	if (m_startDate == value)
		return;
	m_startDate = value;
	emit formChanged();
}

void EphemerisViewModel::setEndDate(const QString& value)
{
	if (m_endDate == value)
		return;
	m_endDate = value;
	emit formChanged();
}

void EphemerisViewModel::setStartTime(const QString& value)
{
	if (m_startTime == value)
		return;
	m_startTime = value;
	emit formChanged();
}

void EphemerisViewModel::setStepDays(double value)
{
	if (m_stepDays == value || (std::isnan(m_stepDays) && std::isnan(value)))
		return;
	m_stepDays = value;
	emit formChanged();
}

void EphemerisViewModel::setPlane(int value)
{
	value = std::clamp(value, 0, 1);
	if (m_plane == value)
		return;
	m_plane = value;
	emit formChanged();
}

void EphemerisViewModel::setScale(int value)
{
	value = std::clamp(value, 0, 2);
	if (m_scale == value)
		return;
	m_scale = value;
	emit formChanged();
}

void EphemerisViewModel::setLightTime(bool value)
{
	if (m_lightTime == value)
		return;
	m_lightTime = value;
	emit formChanged();
}

QString EphemerisViewModel::scaleName() const { return scaleLabel(toScale(m_scale)); }

std::optional<orbital_core::EphemerisRequest> EphemerisViewModel::request(QString* error) const
{
	const auto reject = [error](const QString& message) {
		if (error)
			*error = message;
		return std::optional<orbital_core::EphemerisRequest>{};
	};
	if (!m_solution)
		return reject(tr("No solution"));
	if (!parseTime(m_startTime).isValid())
		return reject(tr("Invalid start time"));
	const auto start = parseMjd(m_startDate, m_startTime);
	if (!start)
		return reject(tr("Invalid start date"));
	const auto end = parseMjd(m_endDate, m_startTime);
	if (!end)
		return reject(tr("Invalid end date"));
	if (!std::isfinite(m_stepDays) || m_stepDays <= 0.0)
		return reject(tr("Step must be positive"));
	if (*end < *start)
		return reject(tr("End is before start"));
	if (m_center == Topocentric && m_station.isEmpty())
		return reject(tr("Enter a station code"));
	std::size_t epochs = 0;
	try {
		epochs = orbital_core::ephemerisEpochCount(*start, *end, m_stepDays);
	} catch (const std::exception&) {
		return reject(tr("Invalid epoch range"));
	}
	if (epochs > orbital_core::kMaxEphemerisEpochs)
		return reject(tr("Too many epochs: %1 (at most %2)").arg(epochs).arg(orbital_core::kMaxEphemerisEpochs));

	orbital_core::EphemerisRequest request;
	request.solution = *m_solution;
	request.center = m_center == Heliocentric ? orbital_core::EphemerisCenter::Heliocentric
	                 : m_center == Geocentric ? orbital_core::EphemerisCenter::Geocentric
	                                          : orbital_core::EphemerisCenter::Station;
	if (request.center == orbital_core::EphemerisCenter::Station)
		request.station = m_station.toStdString();
	request.plane = m_plane == 1 ? orbital_core::EphemerisPlane::Ecliptic : orbital_core::EphemerisPlane::Equatorial;
	request.scale = toScale(m_scale);
	request.startMjd = *start;
	request.endMjd = *end;
	request.stepDays = m_stepDays;
	request.lightTime = m_lightTime;
	return request;
}

int EphemerisViewModel::plannedEpochs() const
{
	const auto start = parseMjd(m_startDate, m_startTime);
	const auto end = parseMjd(m_endDate, m_startTime);
	if (!start || !end || !std::isfinite(m_stepDays) || m_stepDays <= 0.0 || *end < *start)
		return 0;
	try {
		return static_cast<int>(std::min<std::size_t>(orbital_core::ephemerisEpochCount(*start, *end, m_stepDays),
		                                              std::numeric_limits<int>::max()));
	} catch (const std::exception&) {
		return 0;
	}
}

bool EphemerisViewModel::startValid() const { return QDate::fromString(m_startDate.trimmed(), Qt::ISODate).isValid(); }

bool EphemerisViewModel::endValid() const { return QDate::fromString(m_endDate.trimmed(), Qt::ISODate).isValid(); }

bool EphemerisViewModel::timeValid() const { return parseTime(m_startTime).isValid(); }

QString EphemerisViewModel::formError() const
{
	if (!m_solution)
		return {};
	QString error;
	request(&error);
	return error;
}

void EphemerisViewModel::generate()
{
	if (m_busy || m_blocked)
		return;
	QString problem;
	auto task = request(&problem);
	if (!task) {
		m_error = problem;
		emit changed();
		return;
	}
	if (m_worker.joinable())
		m_worker.join();
	const quint64 generation = ++m_generation;
	m_busy = true;
	m_error.clear();
	emit changed();
	m_worker = std::jthread([this, generation, job = std::move(*task)](std::stop_token stop) {
		try {
			auto result =
			    std::make_shared<orbital_core::EphemerisResult>(orbital_core::computeEphemeris(job, {.cancel = stop}));
			QMetaObject::invokeMethod(
			    this, [this, generation, job, result] { finish(generation, job, result); }, Qt::QueuedConnection);
		} catch (const std::exception& e) {
			QMetaObject::invokeMethod(
			    this, [this, generation, message = QString::fromUtf8(e.what())] { fail(generation, message); },
			    Qt::QueuedConnection);
		}
	});
}

void EphemerisViewModel::cancel()
{
	if (m_busy)
		m_worker.request_stop();
}

void EphemerisViewModel::finish(quint64 generation, orbital_core::EphemerisRequest request,
                                std::shared_ptr<orbital_core::EphemerisResult> result)
{
	m_busy = false;
	if (generation == m_generation && !result->cancelled) {
		m_table.setRows(std::move(result->rows), request.plane, request.scale,
		                request.center == orbital_core::EphemerisCenter::Heliocentric, request.solution.magnitudeType);
		result->rows.clear();
		m_lastRequest = std::move(request);
		m_lastResult = std::move(*result);
		m_summary = describe();
	}
	emit changed();
}

void EphemerisViewModel::fail(quint64 generation, const QString& message)
{
	m_busy = false;
	if (generation == m_generation)
		m_error = message;
	emit changed();
}

QString EphemerisViewModel::describe() const
{
	QStringList parts;
	switch (m_lastRequest.center) {
		case orbital_core::EphemerisCenter::Heliocentric:
			parts << tr("Heliocentric");
			break;
		case orbital_core::EphemerisCenter::Geocentric:
			parts << tr("Geocentric");
			break;
		case orbital_core::EphemerisCenter::Station:
			parts << tr("Topocentric %1").arg(QString::fromStdString(m_lastRequest.station));
			break;
	}
	parts << (m_lastRequest.plane == orbital_core::EphemerisPlane::Ecliptic ? tr("ecliptic J2000")
	                                                                        : tr("equatorial J2000"));
	parts << scaleLabel(m_lastRequest.scale);
	if (m_lastRequest.center != orbital_core::EphemerisCenter::Heliocentric)
		parts << (m_lastResult.lightTimeApplied ? tr("light time") : tr("geometric"));
	QString dynamics = m_lastResult.nBody ? tr("N-body") : tr("two-body");
	if (m_lastResult.nonGravitational)
		dynamics += tr(" + non-grav");
	parts << dynamics;
	return parts.join(QStringLiteral(" · "));
}

QString EphemerisViewModel::suggestedFileName(int format) const
{
	static const char* extensions[] = {"txt", "json"};
	QString name = m_designation;
	name.remove(QLatin1Char(' '));
	if (name.isEmpty())
		name = QStringLiteral("object");
	return QStringLiteral("%1_ephemerides.%2").arg(name, QLatin1String(extensions[std::clamp(format, 0, 1)]));
}

bool EphemerisViewModel::save(const QUrl& file, int format)
{
	if (m_busy || !m_table.count() || format < 0 || format > 1)
		return false;
	const QString path = file.isLocalFile() ? file.toLocalFile() : file.toString();
	QSaveFile output(path);
	const QByteArray data = format == 1 ? exportJson() : exportText();
	if (!output.open(QIODevice::WriteOnly) || output.write(data) != data.size() || !output.commit()) {
		m_error = tr("Cannot write %1").arg(path);
		emit changed();
		return false;
	}
	return true;
}

QByteArray EphemerisViewModel::exportText() const
{
	QList<QStringList> lines;
	QStringList header;
	for (int column = 0; column < m_table.columns(); ++column)
		header << m_table.header(column);
	lines << header;
	for (int row = 0; row < m_table.count(); ++row) {
		QStringList cells;
		for (int column = 0; column < m_table.columns(); ++column)
			cells << m_table.cell(row, column);
		lines << cells;
	}
	QList<qsizetype> widths(m_table.columns(), 0);
	for (const auto& line : lines)
		for (int column = 0; column < line.size(); ++column)
			widths[column] = std::max(widths[column], line[column].size());

	QString text;
	text += tr("# Object: %1").arg(m_designation) + QLatin1Char('\n');
	text += tr("# Ephemeris: %1").arg(m_summary) + QLatin1Char('\n');
	text += tr("# Solution epoch: MJD %1 TDB").arg(QString::number(m_lastRequest.solution.epochMjdTdb, 'f', 9)) +
	        QLatin1Char('\n');
	text += tr("# Astrometric positions: light time only, no stellar aberration") + QLatin1Char('\n');
	for (int i = 0; i < lines.size(); ++i) {
		QStringList cells;
		for (int column = 0; column < lines[i].size(); ++column)
			cells << (column == 0 ? lines[i][column].leftJustified(widths[column])
			                      : lines[i][column].rightJustified(widths[column]));
		text += cells.join(QStringLiteral("  ")) + QLatin1Char('\n');
	}
	return text.toUtf8();
}

QByteArray EphemerisViewModel::exportJson() const
{
	QJsonObject root;
	root[QStringLiteral("object")] = m_designation;
	switch (m_lastRequest.center) {
		case orbital_core::EphemerisCenter::Heliocentric:
			root[QStringLiteral("center")] = QStringLiteral("heliocentric");
			break;
		case orbital_core::EphemerisCenter::Geocentric:
			root[QStringLiteral("center")] = QStringLiteral("geocentric");
			break;
		case orbital_core::EphemerisCenter::Station:
			root[QStringLiteral("center")] = QStringLiteral("topocentric");
			root[QStringLiteral("station")] = QString::fromStdString(m_lastRequest.station);
			break;
	}
	root[QStringLiteral("plane")] = m_lastRequest.plane == orbital_core::EphemerisPlane::Ecliptic
	                                    ? QStringLiteral("ecliptic_j2000")
	                                    : QStringLiteral("equatorial_j2000");
	root[QStringLiteral("time_scale")] = scaleLabel(m_lastRequest.scale);
	root[QStringLiteral("light_time")] = m_lastResult.lightTimeApplied;
	root[QStringLiteral("n_body")] = m_lastResult.nBody;
	root[QStringLiteral("non_gravitational")] = m_lastResult.nonGravitational;
	root[QStringLiteral("solution_epoch_mjd_tdb")] = m_lastRequest.solution.epochMjdTdb;
	root[QStringLiteral("absolute_magnitude")] = m_lastRequest.solution.absoluteMagnitude
	                                                 ? QJsonValue(*m_lastRequest.solution.absoluteMagnitude)
	                                                 : QJsonValue(QJsonValue::Null);
	switch (m_lastRequest.solution.magnitudeType) {
		case orbital_core::EphemerisMagnitudeType::CometTotal:
			root[QStringLiteral("magnitude_type")] = QStringLiteral("comet_total");
			root[QStringLiteral("comet_slope_fixed")] = m_lastRequest.solution.cometSlope;
			root[QStringLiteral("comet_phase_coefficient_fixed")] = m_lastRequest.solution.cometPhaseCoefficient;
			break;
		case orbital_core::EphemerisMagnitudeType::CometNuclear:
			root[QStringLiteral("magnitude_type")] = QStringLiteral("comet_nuclear");
			root[QStringLiteral("comet_slope_fixed")] = m_lastRequest.solution.cometSlope;
			root[QStringLiteral("comet_phase_coefficient_fixed")] = m_lastRequest.solution.cometPhaseCoefficient;
			break;
		default:
			root[QStringLiteral("magnitude_type")] = QStringLiteral("asteroid_hg");
			root[QStringLiteral("asteroid_slope_g_fixed")] = m_lastRequest.solution.slopeG;
	}
	QJsonArray rows;
	const QString dateKey = QStringLiteral("date_%1").arg(scaleLabel(m_lastRequest.scale).toLower());
	for (const auto& row : m_table.rows()) {
		QJsonObject item;
		item[dateKey] = EphemerisTableModel::calendar(row.mjd, 3);
		for (const auto& column : numericColumns())
			item[columnKey(column, m_lastRequest.plane)] = jsonNumber(column.value(row));
		rows.append(item);
	}
	root[QStringLiteral("rows")] = rows;
	return QJsonDocument(root).toJson(QJsonDocument::Indented);
}
