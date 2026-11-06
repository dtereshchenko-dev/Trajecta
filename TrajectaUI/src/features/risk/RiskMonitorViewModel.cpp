#include "RiskMonitorViewModel.h"
#include "util/CalculationTime.h"

#include <QDate>
#include <QTime>
#include <QRegularExpression>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QFileInfo>
#include <QUuid>
#include <QStandardPaths>
#include <chrono>
#include <orbital_core/TimeScale.h>
#include <Eigen/Cholesky>
#include <cmath>
#include <limits>

namespace {
	std::optional<double> parseUtc(const QString& value)
	{
		const auto text = value.trimmed();
		static const QRegularExpression pattern(
		    QStringLiteral("^([0-9]{4}-[0-9]{2}-[0-9]{2})[ T]([0-9]{2}:[0-9]{2}:[0-9]{2})(?:\\.([0-9]{1,9}))?$"));
		const auto match = pattern.match(text);
		if (!match.hasMatch())
			return std::nullopt;
		const auto date = QDate::fromString(match.captured(1), Qt::ISODate);
		const auto time = QTime::fromString(match.captured(2), QStringLiteral("HH:mm:ss"));
		if (!date.isValid() || !time.isValid())
			return std::nullopt;
		const double fraction =
		    match.captured(3).isEmpty() ? 0.0 : (QStringLiteral("0.") + match.captured(3)).toDouble();
		return static_cast<double>(date.toJulianDay() - 2400001) +
		       (time.msecsSinceStartOfDay() / 1000.0 + fraction) / 86400.0;
	}

	QString formatUtc(double mjd)
	{
		const double firstDay = static_cast<double>(QDate(1, 1, 1).toJulianDay() - 2400001);
		const double lastDayEnd = static_cast<double>(QDate(9999, 12, 31).toJulianDay() - 2400000);
		if (!std::isfinite(mjd) || mjd < firstDay || mjd >= lastDayEnd)
			return {};
		auto day = static_cast<qint64>(std::floor(mjd));
		auto nanoseconds = static_cast<qint64>(std::llround((mjd - day) * 86400000000000.0));
		if (nanoseconds >= 86400000000000LL) {
			++day;
			nanoseconds = 0;
		}
		const auto date = QDate::fromJulianDay(day + 2400001);
		if (date.year() > 9999)
			return {};
		const auto seconds = nanoseconds / 1000000000LL;
		auto fraction = QString::number(nanoseconds % 1000000000LL).rightJustified(9, QLatin1Char('0'));
		while (fraction.endsWith(QLatin1Char('0')))
			fraction.chop(1);
		return date.toString(QStringLiteral("yyyy-MM-dd")) + QLatin1Char(' ') +
		       QTime(0, 0).addSecs(static_cast<int>(seconds)).toString(QStringLiteral("HH:mm:ss")) +
		       (fraction.isEmpty() ? QString() : QLatin1Char('.') + fraction);
	}

	std::optional<std::uint64_t> parseInteger(const QString& value)
	{
		static const QRegularExpression digits(QStringLiteral("^[0-9]+$"));
		if (!digits.match(value).hasMatch())
			return std::nullopt;
		bool ok = false;
		const auto number = value.toULongLong(&ok);
		return ok ? std::optional<std::uint64_t>(number) : std::nullopt;
	}
}  // namespace

RiskMonitorViewModel::RiskMonitorViewModel(SolutionViewModel& solution, QObject* parent)
    : QObject(parent), m_solution(solution)
{
	connect(&solution, &SolutionViewModel::changed, this, &RiskMonitorViewModel::refreshSource);
	refreshSource();
}

void RiskMonitorViewModel::refreshSource()
{
	const auto designation = m_solution.designation();
	if (!designation.isEmpty()) {
		if (!m_sourceObject.isEmpty() && designation != m_sourceObject) {
			m_centerInitialized = false;
			m_rangeInitialized = false;
			m_centerUtc.clear();
			m_rangeStartUtc.clear();
			m_rangeEndUtc.clear();
		}
		m_sourceObject = designation;
	}
	if (m_worker.joinable())
		m_worker.request_stop();
	m_runs.clear();
	m_error.clear();
	m_selectedRun = -1;
	m_events.setEvents({});
	++m_sourceRevision;
	m_available = false;
	m_covarianceError.clear();
	m_observationError.clear();
	m_observationCount = 0;
	m_radarCount = 0;
	m_parameters = 0;
	const auto& snapshot = m_solution.snapshot();
	if (snapshot && !snapshot->fit.cancelled && snapshot->fit.fitEpochOrbit) {
		const auto& fit = snapshot->fit;
		m_epochMjdTdb = fit.fitEpochOrbit->epochMjd();
		m_available = std::isfinite(m_epochMjdTdb) &&
		              std::isfinite(fit.fitEpochOrbit->getR().dot(fit.fitEpochOrbit->getR())) &&
		              std::isfinite(fit.fitEpochOrbit->getV().dot(fit.fitEpochOrbit->getV()));
		const auto defaultStart = CalculationTime::defaultStartUtc(fit);
		if (m_available && !m_centerInitialized) {
			m_centerUtc = defaultStart ? formatUtc(*defaultStart + m_halfWindowHours / 24.0) : QString();
			m_centerInitialized = !m_centerUtc.isEmpty();
		}
		if (m_available && !m_rangeInitialized) {
			m_rangeStartUtc = defaultStart ? formatUtc(*defaultStart) : QString();
			const double startTdb = defaultStart
			                            ? orbital_core::time::toMjdTdb(*defaultStart, orbital_core::time::Scale::UTC)
			                            : std::numeric_limits<double>::quiet_NaN();
			m_rangeEndUtc =
			    formatUtc(orbital_core::time::fromMjdTdb(startTdb + 30.0 * 365.25, orbital_core::time::Scale::UTC));
			m_rangeInitialized = !m_rangeStartUtc.isEmpty() && !m_rangeEndUtc.isEmpty();
		}
		m_parameters = 6 + fit.finalNonGrav.fitParameterCount();
		const auto& covariance = fit.uncertainty.fitParameterCovariance;
		if (covariance.rows() != m_parameters || covariance.cols() != m_parameters || !covariance.allFinite() ||
		    (covariance.diagonal().array() <= 0.0).any()) {
			m_covarianceError = tr("Full fit-epoch covariance is unavailable.");
		} else {
			const Eigen::VectorXd scales = covariance.diagonal().cwiseSqrt().cwiseInverse();
			const Eigen::MatrixXd correlation = scales.asDiagonal() * covariance * scales.asDiagonal();
			if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-10)) {
				m_covarianceError = tr("Full fit-epoch covariance is not positive definite.");
			} else {
				const Eigen::LLT<Eigen::MatrixXd> factor(correlation);
				if (factor.info() != Eigen::Success)
					m_covarianceError = tr("Full fit-epoch covariance is not positive definite.");
			}
		}
		const auto& all = snapshot->observations.all;
		if (all.empty()) {
			m_observationError = tr("The accepted solution has no observations.");
		} else {
			bool indicesValid = true;
			std::vector<bool> seen(all.size());
			for (const auto index : fit.keptObsIndices) {
				if (index >= all.size() || seen[index]) {
					indicesValid = false;
					break;
				}
				seen[index] = true;
			}
			if (!indicesValid) {
				m_observationError = tr("The accepted observation mask is invalid.");
			} else {
				for (const auto& observation : snapshot->fittedObservations()) {
					++m_observationCount;
					if (observation.radar)
						++m_radarCount;
				}
			}
		}
	}
	emit changed();
}

void RiskMonitorViewModel::setMethod(int value)
{
	if (value < 0 || value > 1 || value == m_method)
		return;
	m_method = value;
	emit changed();
}

void RiskMonitorViewModel::setTarget(int value)
{
	if (value < 0 || value > 2 || value == m_target)
		return;
	m_target = value;
	emit changed();
}

void RiskMonitorViewModel::setIntervalMode(int value)
{
	if (value < 0 || value > 1 || value == m_intervalMode)
		return;
	m_intervalMode = value;
	emit changed();
}

void RiskMonitorViewModel::setCenterUtc(const QString& value)
{
	m_centerInitialized = true;
	if (value == m_centerUtc)
		return;
	m_centerUtc = value;
	emit changed();
}

void RiskMonitorViewModel::setHalfWindowHours(double value)
{
	if (value == m_halfWindowHours || (std::isnan(value) && std::isnan(m_halfWindowHours)))
		return;
	m_halfWindowHours = value;
	emit changed();
}

void RiskMonitorViewModel::setRangeStartUtc(const QString& value)
{
	if (value == m_rangeStartUtc)
		return;
	m_rangeInitialized = true;
	m_rangeStartUtc = value;
	emit changed();
}

void RiskMonitorViewModel::setRangeEndUtc(const QString& value)
{
	if (value == m_rangeEndUtc)
		return;
	m_rangeInitialized = true;
	m_rangeEndUtc = value;
	emit changed();
}

double RiskMonitorViewModel::rangeYears() const
{
	const auto start = parseUtc(m_rangeStartUtc);
	const auto end = parseUtc(m_rangeEndUtc);
	if (!start || !end)
		return std::numeric_limits<double>::quiet_NaN();
	const double startTdb = orbital_core::time::toMjdTdb(*start, orbital_core::time::Scale::UTC);
	const double endTdb = orbital_core::time::toMjdTdb(*end, orbital_core::time::Scale::UTC);
	return (endTdb - startTdb) / 365.25;
}

void RiskMonitorViewModel::setRangeYears(double value)
{
	const auto start = parseUtc(m_rangeStartUtc);
	if (!start || !std::isfinite(value) || value <= 0.0) {
		setRangeEndUtc({});
		return;
	}
	const double startTdb = orbital_core::time::toMjdTdb(*start, orbital_core::time::Scale::UTC);
	setRangeEndUtc(
	    formatUtc(orbital_core::time::fromMjdTdb(startTdb + value * 365.25, orbital_core::time::Scale::UTC)));
}

void RiskMonitorViewModel::setSamples(const QString& value)
{
	if (value == m_samples[m_method])
		return;
	m_samples[m_method] = value;
	emit changed();
}

void RiskMonitorViewModel::setSeed(const QString& value)
{
	if (value == m_seed)
		return;
	m_seed = value;
	emit changed();
}

void RiskMonitorViewModel::setBoundary(int value)
{
	if (value < 0 || value > 1 || value == m_boundary || m_target != 0)
		return;
	m_boundary = value;
	emit changed();
}

QString RiskMonitorViewModel::sourceError() const
{
	if (!m_available)
		return tr("A nominal state at the fit epoch is required.");
	if (!m_solution.snapshot()->fit.converged)
		return tr("A converged solution is required.");
	if (!m_solution.snapshot()->effectiveConfig.nbody.enabled)
		return tr("An N-body solution is required.");
	return m_method == 0 ? m_covarianceError : m_observationError;
}

QString RiskMonitorViewModel::formError() const
{
	if (m_intervalMode == 0) {
		if (!parseUtc(m_centerUtc))
			return tr("Enter the window centre as YYYY-MM-DD HH:MM:SS UTC.");
		if (!std::isfinite(m_halfWindowHours) || m_halfWindowHours <= 0.0)
			return tr("The half-window must be a finite positive number of hours.");
	} else {
		if (!parseUtc(m_rangeStartUtc))
			return tr("Enter the range start as YYYY-MM-DD HH:MM:SS UTC.");
		if (!parseUtc(m_rangeEndUtc))
			return tr("Enter the range end as YYYY-MM-DD HH:MM:SS UTC.");
	}
	const auto [start, end] = *intervalBounds();
	if (end <= start)
		return tr("The calculation interval must end after its start.");
	const double firstDay = static_cast<double>(QDate(1, 1, 1).toJulianDay() - 2400001);
	const double lastDayEnd = static_cast<double>(QDate(9999, 12, 31).toJulianDay() - 2400000);
	if (!std::isfinite(start) || !std::isfinite(end) || start < firstDay || end >= lastDayEnd)
		return tr("The calculation interval is outside the supported calendar range.");
	if (m_available) {
		try {
			const double startTdb = orbital_core::time::toMjdTdb(start, orbital_core::time::Scale::UTC);
			if (!std::isfinite(startTdb))
				return tr("UTC conversion is unavailable for this window.");
			if (startTdb < m_epochMjdTdb)
				return tr("The calculation interval must start at or after the fit epoch.");
			if (!std::isfinite(orbital_core::time::toMjdTdb(end, orbital_core::time::Scale::UTC)))
				return tr("UTC conversion is unavailable for this window.");
		} catch (...) {
			return tr("UTC conversion is unavailable for this window.");
		}
	}
	const auto count = parseInteger(samples());
	if (!count || *count == 0 || *count > 1000000)
		return tr("The sample count must be an integer from 1 to 1000000.");
	if (!parseInteger(m_seed))
		return tr("Seed must be an integer from 0 to 18446744073709551615.");
	return {};
}

QString RiskMonitorViewModel::fitEpoch() const
{
	return m_available ? QString::number(m_epochMjdTdb, 'f', 9) + QStringLiteral(" TDB") : QString();
}

QString RiskMonitorViewModel::covarianceSize() const
{
	return m_available && m_covarianceError.isEmpty() ? QStringLiteral("%1 x %1").arg(m_parameters) : tr("Unavailable");
}

QString RiskMonitorViewModel::windowDuration() const
{
	if (m_intervalMode == 1) {
		const double years = rangeYears();
		return std::isfinite(years) && years > 0.0 ? tr("%1 years").arg(QString::number(years, 'g', 8)) : QString();
	}
	if (!std::isfinite(m_halfWindowHours) || m_halfWindowHours <= 0.0)
		return {};
	return tr("%1 days").arg(QString::number(m_halfWindowHours / 12.0, 'g', 8));
}

std::optional<std::pair<double, double>> RiskMonitorViewModel::intervalBounds() const
{
	if (m_intervalMode == 1) {
		const auto start = parseUtc(m_rangeStartUtc);
		const auto end = parseUtc(m_rangeEndUtc);
		if (!start || !end)
			return std::nullopt;
		return std::pair{*start, *end};
	}
	const auto center = parseUtc(m_centerUtc);
	if (!center)
		return std::nullopt;
	return std::pair{*center - m_halfWindowHours / 24.0, *center + m_halfWindowHours / 24.0};
}

std::optional<RiskMonitorDraft> RiskMonitorViewModel::draft() const
{
	if (!inputReady())
		return std::nullopt;
	const auto [start, end] = *intervalBounds();
	return RiskMonitorDraft{static_cast<RiskMonitorDraft::Method>(m_method),
	                        static_cast<RiskMonitorDraft::Target>(m_target),
	                        static_cast<RiskMonitorDraft::IntervalMode>(m_intervalMode),
	                        start,
	                        end,
	                        altitudeKm(),
	                        static_cast<int>(*parseInteger(samples())),
	                        *parseInteger(m_seed),
	                        m_sourceRevision};
}

RiskMonitorViewModel::~RiskMonitorViewModel()
{
	if (m_worker.joinable()) {
		m_worker.request_stop();
		m_worker.join();
	}
}

void RiskMonitorViewModel::setBlocked(bool value)
{
	if (m_blocked == value)
		return;
	m_blocked = value;
	emit changed();
}

QStringList RiskMonitorViewModel::runs() const
{
	QStringList result;
	for (const auto& run : m_runs)
		result.push_back(run.label);
	return result;
}

void RiskMonitorViewModel::setSelectedRun(int value)
{
	if (value < 0 || value >= static_cast<int>(m_runs.size()) || value == m_selectedRun)
		return;
	m_selectedRun = value;
	m_events.setEvents(m_runs[value].result->events);
	emit changed();
}

QVariantList RiskMonitorViewModel::atmosphericSources() const
{
	QVariantList sources;
	for (const auto& run : m_runs) {
		const auto& r = *run.result;
		if (r.config.target != orbital_core::hazard::Planet::Earth || r.config.altitudeKm != 100.0 || r.cancelled ||
		    r.completed != r.config.samples || r.entryClones == 0 || run.archive.isEmpty())
			continue;
		sources.push_back(QVariantMap{{"label", run.label},
		                              {"path", run.archive + "/result.json"},
		                              {"url", QUrl::fromLocalFile(run.archive + "/result.json")},
		                              {"entries", r.entryClones}});
	}
	return sources;
}

QString RiskMonitorViewModel::summary() const
{
	if (m_selectedRun < 0)
		return {};
	const auto& run = m_runs[m_selectedRun];
	const auto& r = *run.result;
	return tr("Complete %1 / %2  |  Fit errors %3  |  Propagation errors %4  |  Stopped %5  |  Not started %6  |  %7 s")
	    .arg(r.completed)
	    .arg(r.config.samples)
	    .arg(r.failedFit)
	    .arg(r.failedPropagation)
	    .arg(r.cancelledClones)
	    .arg(r.unstarted)
	    .arg(run.seconds, 0, 'f', 1);
}

QString RiskMonitorViewModel::runDetails() const
{
	if (m_selectedRun < 0)
		return {};
	const auto& config = m_runs[m_selectedRun].result->config;
	const auto date = [](double tdb) {
		return formatUtc(orbital_core::time::fromMjdTdb(tdb, orbital_core::time::Scale::UTC));
	};
	return tr("UTC %1 to %2  |  Altitude %3 km  |  Seed %4  |  RA15")
	    .arg(date(config.startMjdTdb), date(config.endMjdTdb))
	    .arg(config.altitudeKm, 0, 'g', 8)
	    .arg(config.seed);
}

QString RiskMonitorViewModel::warning() const
{
	if (m_selectedRun < 0)
		return {};
	const auto& r = *m_runs[m_selectedRun].result;
	QStringList parts;
	if (r.cancelled)
		parts << tr("Cancelled: partial result.");
	if (r.completed != r.config.samples)
		parts << tr("Incomplete sample: probabilities are conditional on successful trials and may be biased.");
	if (r.nominal.status != orbital_core::hazard::McCloneStatus::Complete)
		parts << tr("Nominal propagation failed: %1").arg(QString::fromStdString(r.nominal.error));
	if (r.entriesBeforeInterval)
		parts << tr("%1 clones crossed the boundary before the requested interval.").arg(r.entriesBeforeInterval);
	return parts.join(QLatin1Char(' '));
}

QString RiskMonitorViewModel::archivePath() const
{
	return m_selectedRun < 0 ? QString() : m_runs[m_selectedRun].archive;
}

void RiskMonitorViewModel::openArchive(bool report)
{
	auto path = archivePath();
	if (path.isEmpty())
		return;
	if (report)
		path = QDir(path).filePath(QStringLiteral("report.html"));
	if (!QDesktopServices::openUrl(QUrl::fromLocalFile(path))) {
		m_error = tr("Cannot open the MC archive.");
		emit changed();
	}
}

void RiskMonitorViewModel::cancel()
{
	if (!m_busy || m_cancelling)
		return;
	m_cancelling = true;
	m_progress = tr("Cancelling...");
	m_worker.request_stop();
	emit changed();
}

void RiskMonitorViewModel::run()
{
	if (!canRun())
		return;
	const auto request = *draft();
	auto source = *m_solution.snapshot();
	source.fit.fitEpochOrbit = std::make_shared<Orbit>(*source.fit.fitEpochOrbit);
	if (source.fit.orbit)
		source.fit.orbit = std::make_shared<Orbit>(*source.fit.orbit);
	const QString root =
	    QDir(m_solution.directory().isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)
	                                          : m_solution.directory())
	        .filePath(QStringLiteral("MC"));
	const QString label = QDateTime::currentDateTimeUtc().toString(QStringLiteral("yyyy-MM-dd HH:mm:ss 'UTC'")) +
	                      QStringLiteral("  |  ") + (m_method == 0 ? tr("Covariance MC") : tr("Observation MC")) +
	                      QStringLiteral("  |  ") + QStringList{tr("Earth"), tr("Moon"), tr("Mars")}[m_target] +
	                      QStringLiteral("  |  N=") + QString::number(request.samples);
	m_busy = true;
	m_cancelling = false;
	m_progressFraction = 0;
	m_progress = tr("Preparing MC...");
	m_error.clear();
	emit changed();
	emit runStarted();
	m_worker = std::jthread([this, source = std::move(source), request, root, label](std::stop_token stop) {
		namespace hz = orbital_core::hazard;
		const auto started = std::chrono::steady_clock::now();
		std::shared_ptr<hz::IntervalMcResult> result;
		QString archive, error, stagingArchive;
		try {
			hz::IntervalMcConfig config;
			config.method = request.method == RiskMonitorDraft::Method::Covariance ? hz::McMethod::Covariance
			                                                                       : hz::McMethod::Observations;
			config.target = request.target == RiskMonitorDraft::Target::Earth  ? hz::Planet::Earth
			                : request.target == RiskMonitorDraft::Target::Moon ? hz::Planet::Moon
			                                                                   : hz::Planet::Mars;
			config.startMjdTdb = orbital_core::time::toMjdTdb(request.startMjdUtc, orbital_core::time::Scale::UTC);
			config.endMjdTdb = orbital_core::time::toMjdTdb(request.endMjdUtc, orbital_core::time::Scale::UTC);
			config.altitudeKm = request.altitudeKm;
			config.samples = request.samples;
			config.seed = request.seed;
			auto lastProgress = started;
			result = std::make_shared<hz::IntervalMcResult>(
			    hz::runIntervalMonteCarlo(source, config, stop, [this, &lastProgress](int done, int total) {
				    const auto now = std::chrono::steady_clock::now();
				    if (done != total && now - lastProgress < std::chrono::milliseconds(100))
					    return;
				    lastProgress = now;
				    QMetaObject::invokeMethod(
				        this,
				        [this, done, total] {
					        if (!m_cancelling)
						        m_progress = tr("Monte Carlo: %1 / %2").arg(done).arg(total);
					        m_progressFraction = double(done) / total;
					        emit changed();
				        },
				        Qt::QueuedConnection);
			    }));
			if (!QDir().mkpath(root))
				throw std::runtime_error("Cannot create the MC archive folder");
			const auto path =
			    QDir(root).filePath(QStringLiteral("run-") + QUuid::createUuid().toString(QUuid::WithoutBraces));
			stagingArchive = path + QStringLiteral(".partial");
			hz::saveIntervalMonteCarlo(std::filesystem::path(path.toStdWString()), source, *result);
			archive = path;
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
			if (!stagingArchive.isEmpty() &&
			    QFileInfo::exists(QDir(stagingArchive).filePath(QStringLiteral("result.json"))) &&
			    QFileInfo::exists(QDir(stagingArchive).filePath(QStringLiteral("report.html"))))
				archive = stagingArchive;
		}
		if (result) {
			result->clones.clear();
			result->clones.shrink_to_fit();
		}
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
		QMetaObject::invokeMethod(
		    this,
		    [this, result, archive, error, label, seconds, request] {
			    m_busy = false;
			    m_cancelling = false;
			    if (request.sourceRevision == m_sourceRevision) {
				    m_error = error;
				    if (result) {
					    m_runs.push_back({result, label, archive, seconds});
					    m_selectedRun = -1;
					    setSelectedRun(static_cast<int>(m_runs.size()) - 1);
					    m_progress = result->cancelled ? tr("MC cancelled") : tr("MC completed");
				    } else
					    m_progress = tr("MC failed");
			    }
			    emit changed();
		    },
		    Qt::QueuedConnection);
	});
}
