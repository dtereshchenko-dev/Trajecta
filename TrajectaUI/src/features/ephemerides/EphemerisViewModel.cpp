#include "EphemerisViewModel.h"
#include "util/CalculationTime.h"

#include <QDate>
#include <QDesktopServices>
#include <QDir>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QTime>
#include <QUuid>
#include <orbital_core/Types.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace {

	constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;
	constexpr double kRadPerDayToArcsecPerHour = 206264.80624709636 / 24.0;
	constexpr double kAuPerDayToKmPerSecond = 149597870.7 / 86400.0;

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

	QString textColumn(const QString& value)
	{
		QString text = value;
		return text.replace(QStringLiteral("λ"), QStringLiteral("lon"))
		    .replace(QStringLiteral("β"), QStringLiteral("lat"))
		    .replace(QStringLiteral("Δ"), QStringLiteral("Delta"))
		    .replace(QStringLiteral("σ"), QStringLiteral("sigma"))
		    .replace(QStringLiteral("²"), QStringLiteral("^2"))
		    .replace(QStringLiteral("°"), QStringLiteral("deg"))
		    .replace(QStringLiteral("′"), QStringLiteral("'"))
		    .replace(QStringLiteral("″"), QStringLiteral("\""))
		    .replace(QStringLiteral("—"), QStringLiteral("-"));
	}

	QJsonValue jsonNumber(double value)
	{
		return std::isfinite(value) ? QJsonValue(value) : QJsonValue(QJsonValue::Null);
	}

	struct NumericColumn {
		const char* key;
		EphemerisTableModel::Column column;
		double (*value)(const orbital_core::EphemerisRow&);
	};

	const std::vector<NumericColumn>& numericColumns()
	{
		using Row = orbital_core::EphemerisRow;
		using C = EphemerisTableModel;
		static const std::vector<NumericColumn> columns{
		    {"mjd_utc", C::Date, [](const Row& r) { return r.mjdUtc; }},
		    {"mjd_tt", C::Date, [](const Row& r) { return r.mjdTt; }},
		    {"mjd_tdb", C::Date, [](const Row& r) { return r.mjdTdb; }},
		    {"lon_deg", C::Longitude, [](const Row& r) { return r.longitude * kRadToDeg; }},
		    {"lat_deg", C::Latitude, [](const Row& r) { return r.latitude * kRadToDeg; }},
		    {"lon_rate_coslat_arcsec_h", C::LongitudeRate,
		     [](const Row& r) { return r.longitudeRate * kRadPerDayToArcsecPerHour; }},
		    {"lat_rate_arcsec_h", C::LatitudeRate,
		     [](const Row& r) { return r.latitudeRate * kRadPerDayToArcsecPerHour; }},
		    {"delta_au", C::Distance, [](const Row& r) { return r.distanceAu; }},
		    {"r_au", C::SolarDistance, [](const Row& r) { return r.heliocentricAu; }},
		    {"delta_rate_km_s", C::DistanceRate,
		     [](const Row& r) { return r.distanceRateAuPerDay * kAuPerDayToKmPerSecond; }},
		    {"r_rate_km_s", C::SolarDistanceRate,
		     [](const Row& r) { return r.heliocentricRateAuPerDay * kAuPerDayToKmPerSecond; }},
		    {"one_way_light_time_min", C::LightTime, [](const Row& r) { return r.oneWayLightTimeDays * 1440.0; }},
		    {"elongation_deg", C::Elongation, [](const Row& r) { return r.elongation * kRadToDeg; }},
		    {"phase_deg", C::Phase, [](const Row& r) { return r.phase * kRadToDeg; }},
		    {"sky_motion_arcsec_h", C::SkyMotion, &C::skyMotion},
		    {"sky_motion_pa_deg", C::SkyPositionAngle, &C::skyPositionAngle},
		    {"lon_coslat_3sigma_arcsec", C::LongitudeUncertainty,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->longitude3SigmaArcsec : std::numeric_limits<double>::quiet_NaN(); }},
		    {"lat_3sigma_arcsec", C::LatitudeUncertainty,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->latitude3SigmaArcsec : std::numeric_limits<double>::quiet_NaN(); }},
		    {"ellipse_major_3sigma_arcsec", C::EllipseMajor,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->major3SigmaArcsec : std::numeric_limits<double>::quiet_NaN(); }},
		    {"ellipse_minor_3sigma_arcsec", C::EllipseMinor,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->minor3SigmaArcsec : std::numeric_limits<double>::quiet_NaN(); }},
		    {"ellipse_angle_deg", C::EllipseAngle,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->positionAngleDeg : std::numeric_limits<double>::quiet_NaN(); }},
		    {"ellipse_area_3sigma_arcsec2", C::EllipseArea,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->area3SigmaArcsec2 : std::numeric_limits<double>::quiet_NaN(); }},
		    {"sky_3sigma_arcsec", C::SkyUncertainty,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->sky3SigmaArcsec : std::numeric_limits<double>::quiet_NaN(); }},
		    {"range_3sigma_km", C::RangeUncertainty,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->range3SigmaKm : std::numeric_limits<double>::quiet_NaN(); }},
		    {"range_rate_3sigma_km_s", C::RangeRateUncertainty,
		     [](const Row& r) { return r.uncertainty ? r.uncertainty->rangeRate3SigmaKmPerSecond : std::numeric_limits<double>::quiet_NaN(); }},
		    {"v_mag", C::Magnitude,
		     [](const Row& r) { return r.visualMagnitude.value_or(std::numeric_limits<double>::quiet_NaN()); }},
		    {"x_au", C::X, [](const Row& r) { return r.position.x; }},
		    {"y_au", C::Y, [](const Row& r) { return r.position.y; }},
		    {"z_au", C::Z, [](const Row& r) { return r.position.z; }},
		    {"vx_au_d", C::Vx, [](const Row& r) { return r.velocity.x; }},
		    {"vy_au_d", C::Vy, [](const Row& r) { return r.velocity.y; }},
		    {"vz_au_d", C::Vz, [](const Row& r) { return r.velocity.z; }},
		};
		return columns;
	}

	QString columnKey(const NumericColumn& column, orbital_core::EphemerisPlane plane)
	{
		QString key = QString::fromLatin1(column.key);
		if (plane == orbital_core::EphemerisPlane::Equatorial)
			key.replace(QStringLiteral("lon_rate_coslat"), QStringLiteral("ra_rate_cosdec"))
			    .replace(QStringLiteral("lon_coslat_3sigma"), QStringLiteral("ra_cosdec_3sigma"))
			    .replace(QStringLiteral("lat_3sigma"), QStringLiteral("dec_3sigma"))
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
                                     const orbital_core::photometry::CometSolution* cometPhotometry,
                                     const QString& directory)
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
	const auto& uncertainty = result.fit.uncertainty;
	const int dimension = 6 + solution.nonGrav.fitParameterCount();
	if (result.fit.fitEpochOrbit && uncertainty.fitParameterCovariance.rows() == dimension &&
	    uncertainty.fitParameterCovariance.cols() == dimension)
		solution.covariance = uncertainty.fitParameterCovariance;
	else if (result.fit.fitEpochOrbit && dimension == 6 && uncertainty.hasStateCovariance())
		solution.covariance = uncertainty.cartesianStateCovariance;
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
	m_solutionConverged = result.fit.converged;

	if (designation != m_designation || m_startDate.isEmpty()) {
		const auto startUtc = CalculationTime::defaultStartUtc(result.fit);
		m_startDate = startUtc ? isoDate(*startUtc) : QString();
		m_endDate = startUtc ? isoDate(std::floor(*startUtc) + 14.0) : QString();
	}
	m_designation = designation;
	m_directory = directory;
	m_archivePath.clear();
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
	m_archivePath.clear();
	m_directory.clear();
	m_solutionConverged = false;
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
	if (!m_solutionConverged)
		return reject(tr("A converged solution is required."));
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
	request.uncertainty = true;
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
	const QString designation = m_designation;
	const QString root = QDir(m_directory.isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)
	                                                : m_directory).filePath(QStringLiteral("Ephemerides"));
	const auto layout = m_table.layout();
	std::vector<std::pair<int, bool>> columns;
	for (const auto& option : m_table.columnOptions()) {
		const auto values = option.toMap();
		columns.emplace_back(values.value("id").toInt(), values.value("enabled").toBool());
	}
	m_busy = true;
	m_error.clear();
	emit changed();
	m_worker = std::jthread([this, generation, job = std::move(*task), designation, root, layout, columns](std::stop_token stop) {
		try {
			auto result =
			    std::make_shared<orbital_core::EphemerisResult>(orbital_core::computeEphemeris(job, {.cancel = stop}));
			QString archivePath, archiveError;
			if (!result->cancelled && !stop.stop_requested()) {
				EphemerisViewModel archive;
				archive.m_designation = designation;
				archive.m_lastRequest = job;
				archive.m_table.setLayout(layout);
				for (const auto& [column, enabled] : columns)
					archive.m_table.setColumnEnabled(column, enabled);
				archive.m_table.setRows(std::move(result->rows), job.plane, job.scale,
				                        job.center == orbital_core::EphemerisCenter::Heliocentric, job.solution.magnitudeType);
				archive.m_lastResult = *result;
				archive.m_summary = archive.describe();
				try {
					if (!QDir().mkpath(root))
						throw std::runtime_error("Cannot create ephemeris archive folder");
					QTemporaryDir staging(QDir(root).filePath(QStringLiteral(".staging-XXXXXX")));
					if (!staging.isValid())
						throw std::runtime_error("Cannot create temporary ephemeris archive");
					if (!archive.save(QUrl::fromLocalFile(staging.filePath(QStringLiteral("ephemeris.txt"))), 0) ||
					    !archive.save(QUrl::fromLocalFile(staging.filePath(QStringLiteral("result.json"))), 1))
						throw std::runtime_error(archive.error().toStdString());
					if (stop.stop_requested()) {
						result->cancelled = true;
					} else {
						const QString path = QDir(root).filePath(QStringLiteral("run-") + QUuid::createUuid().toString(QUuid::WithoutBraces));
						if (!QDir().rename(staging.path(), path))
							throw std::runtime_error("Cannot finalize ephemeris archive");
						staging.setAutoRemove(false);
						archivePath = path;
					}
				} catch (const std::exception& e) {
					archiveError = QString::fromUtf8(e.what());
				}
				result->rows = archive.m_table.takeRows();
			} else if (stop.stop_requested()) {
				result->cancelled = true;
			}
			QMetaObject::invokeMethod(
			    this, [this, generation, job, result, archivePath, archiveError] {
				    finish(generation, job, result, archivePath, archiveError);
			    }, Qt::QueuedConnection);
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
	                            std::shared_ptr<orbital_core::EphemerisResult> result, const QString& archive,
	                            const QString& archiveError)
{
	m_busy = false;
	if (generation == m_generation && !result->cancelled) {
		m_table.setRows(std::move(result->rows), request.plane, request.scale,
		                request.center == orbital_core::EphemerisCenter::Heliocentric, request.solution.magnitudeType);
		result->rows.clear();
		m_lastRequest = std::move(request);
		m_lastResult = std::move(*result);
		m_summary = describe();
		m_archivePath = archive;
		if (!archiveError.isEmpty())
			m_error = tr("Ephemeris calculated, but could not be saved: %1").arg(archiveError);
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
			parts << (m_lastRequest.station == "500" ? tr("Geocentric (MPC 500)")
			                                        : tr("Topocentric %1").arg(QString::fromStdString(m_lastRequest.station)));
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

bool EphemerisViewModel::openArchive()
{
	if (!m_archivePath.isEmpty() && QDesktopServices::openUrl(QUrl::fromLocalFile(m_archivePath)))
		return true;
	m_error = tr("Cannot open the ephemeris archive folder.");
	emit changed();
	return false;
}

bool EphemerisViewModel::save(const QUrl& file, int format)
{
	if (m_busy || !m_table.count() || format < 0 || format > 1)
		return false;
	if (!file.isLocalFile()) {
		m_error = tr("Choose a local file to save ephemerides.");
		emit changed();
		return false;
	}
	const QString path = file.toLocalFile();
	QSaveFile output(path);
	const QByteArray data = format == 1 ? exportJson() : exportText();
	if (!output.open(QIODevice::WriteOnly) || output.write(data) != data.size() || !output.commit()) {
		m_error = tr("Cannot write %1: %2").arg(QDir::toNativeSeparators(path), output.errorString());
		emit changed();
		return false;
	}
	if (!m_error.isEmpty()) {
		m_error.clear();
		emit changed();
	}
	return true;
}

QByteArray EphemerisViewModel::exportText() const
{
	using Column = EphemerisTableModel::Column;
	QList<Column> columns;
	for (int column = 0; column < m_table.columns(); ++column) {
		const int logical = m_table.logicalColumn(column);
		if (logical >= Column::Date && logical < Column::ColumnCount)
			columns << static_cast<Column>(logical);
	}
	for (int column = Column::LongitudeUncertainty; column <= Column::RangeRateUncertainty; ++column)
		if (!columns.contains(static_cast<Column>(column)))
			columns << static_cast<Column>(column);
	QList<QStringList> lines;
	QStringList header;
	for (const auto column : columns)
		header << textColumn(m_table.headerForColumn(column));
	lines << header;
	for (int row = 0; row < m_table.count(); ++row) {
		QStringList cells;
		for (const auto column : columns)
			cells << textColumn(m_table.cellForColumn(row, column));
		lines << cells;
	}
	QList<qsizetype> widths(columns.size(), 0);
	for (const auto& line : lines)
		for (int column = 0; column < line.size(); ++column)
			widths[column] = std::max(widths[column], line[column].size());

	QString text;
	text += tr("# Object: %1").arg(m_designation) + QLatin1Char('\n');
	text += tr("# Ephemeris: %1").arg(m_summary) + QLatin1Char('\n');
	text += tr("# Solution epoch: MJD %1 TDB").arg(QString::number(m_lastRequest.solution.epochMjdTdb, 'f', 9)) +
	        QLatin1Char('\n');
	text += m_lastResult.lightTimeApplied
	            ? tr("# Positions: light time applied; no stellar aberration") + QLatin1Char('\n')
	            : tr("# Positions: geometric; no light time or stellar aberration") + QLatin1Char('\n');
	if (m_table.includesColumn(EphemerisTableModel::Elongation) &&
	    m_lastRequest.center != orbital_core::EphemerisCenter::Heliocentric)
		text += m_lastResult.lightTimeApplied
		            ? tr("# Elongation: apparent directions with light time, solar deflection and stellar aberration.") +
		                  QLatin1Char('\n')
		            : tr("# Elongation: geometric Sun and target directions.") + QLatin1Char('\n');
	text += tr("# Light time: displayed range divided by the speed of light.") + QLatin1Char('\n');
	text += tr("# Range rates: radial relative velocities in km/s.") + QLatin1Char('\n');
	if (m_lastResult.lightTimeApplied) {
		text += tr("# Range rates are not derivatives of light-time-corrected ranges.") + QLatin1Char('\n');
		text +=
		    tr("# Solar range includes Sun-to-target light time; r rate is at target emission.") + QLatin1Char('\n');
	}
	if (m_table.includesColumn(EphemerisTableModel::SkyPositionAngle))
		text += tr("# Sky motion PA is measured from north towards east in the selected plane.") + QLatin1Char('\n');
	text += tr("# Uncertainties: linearized 3-sigma fit covariance propagated from the fit epoch.") + QLatin1Char('\n');
	text += tr("# Ellipse angle: from increasing longitude toward increasing latitude, modulo 180 deg.") + QLatin1Char('\n');
	if (!m_lastResult.uncertaintyComputed)
		text += tr("# Uncertainties unavailable: fit covariance missing, invalid, or not computed.") + QLatin1Char('\n');
	for (int i = 0; i < lines.size(); ++i) {
		if (i == 1) {
			QStringList separators;
			for (const auto width : widths)
				separators << QString(width, QLatin1Char('-'));
			text += separators.join(QStringLiteral("   ")) + QLatin1Char('\n');
			text += QStringLiteral("$$SOE\n");
		}
		QStringList cells;
		for (int column = 0; column < lines[i].size(); ++column)
			cells << (column == 0 ? lines[i][column].leftJustified(widths[column])
			                      : lines[i][column].rightJustified(widths[column]));
		text += cells.join(QStringLiteral("   ")) + QLatin1Char('\n');
	}
	text += QStringLiteral("$$EOE\n");
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
			root[QStringLiteral("center")] = m_lastRequest.station == "500" ? QStringLiteral("geocentric")
			                                                                : QStringLiteral("topocentric");
			root[QStringLiteral("station")] = QString::fromStdString(m_lastRequest.station);
			break;
	}
	root[QStringLiteral("plane")] = m_lastRequest.plane == orbital_core::EphemerisPlane::Ecliptic
	                                    ? QStringLiteral("ecliptic_j2000")
	                                    : QStringLiteral("equatorial_j2000");
	root[QStringLiteral("layout")] = m_table.layout() == EphemerisTableModel::Layout::Cartesian
	                                     ? QStringLiteral("cartesian")
	                                     : QStringLiteral("spherical");
	if (m_table.includesColumn(EphemerisTableModel::SkyPositionAngle))
		root[QStringLiteral("sky_motion_pa_convention")] = QStringLiteral("north_through_east");
	root[QStringLiteral("time_scale")] = scaleLabel(m_lastRequest.scale);
	root[QStringLiteral("light_time")] = m_lastResult.lightTimeApplied;
	root[QStringLiteral("uncertainty_definition")] = QStringLiteral("linearized_3sigma_fit_covariance");
	root[QStringLiteral("uncertainty_computed")] = m_lastResult.uncertaintyComputed;
	root[QStringLiteral("ellipse_angle_convention")] = QStringLiteral("longitude_toward_latitude_modulo_180_deg");
	if (m_table.includesColumn(EphemerisTableModel::Elongation) &&
	    m_lastRequest.center != orbital_core::EphemerisCenter::Heliocentric)
		root[QStringLiteral("elongation_definition")] =
		    m_lastResult.lightTimeApplied ? QStringLiteral("apparent") : QStringLiteral("geometric");
	root[QStringLiteral("range_rate_definition")] = QStringLiteral("radial_relative_velocity");
	root[QStringLiteral("solar_range_light_time")] = m_lastResult.lightTimeApplied;
	root[QStringLiteral("light_time_definition")] = QStringLiteral("displayed_range_divided_by_c");
	root[QStringLiteral("n_body")] = m_lastResult.nBody;
	const auto& forces = m_lastRequest.solution.forces;
	root[QStringLiteral("forces")] = QJsonObject{
	    {"planet_mask", forces.planetMask}, {"sb441", forces.useSb441},
	    {"relativity", orbital_core::relativityModeName(forces.forceModel.relativity)},
	    {"sun_j2", true}, {"earth_j2", forces.forceModel.earthJ2},
	    {"earth_j3", forces.forceModel.earthJ3}, {"earth_j4", forces.forceModel.earthJ4},
	    {"earth_rotating_pole", forces.forceModel.rotatingEarthPole}};
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
			if (m_table.includesColumn(column.column) ||
			    (column.column >= EphemerisTableModel::LongitudeUncertainty &&
			     column.column <= EphemerisTableModel::RangeRateUncertainty))
				item[columnKey(column, m_lastRequest.plane)] = jsonNumber(column.value(row));
		rows.append(item);
	}
	root[QStringLiteral("rows")] = rows;
	return QJsonDocument(root).toJson(QJsonDocument::Indented);
}
