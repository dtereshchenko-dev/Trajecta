#include "CloseApproachesViewModel.h"
#include <orbital_core/TimeScale.h>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QRegularExpression>
#include <QSaveFile>
#include <QStandardPaths>
#include <QUuid>
#include <chrono>
#include <cmath>

namespace {
	namespace hz = orbital_core::hazard;

	std::optional<double> utcToTdb(const QString& value)
	{
		static const QRegularExpression pattern(
		    QStringLiteral("^([0-9]{4}-[0-9]{2}-[0-9]{2})[ T]([0-9]{2}:[0-9]{2}:[0-9]{2})(?:\\.([0-9]{1,9}))?$"));
		const auto match = pattern.match(value.trimmed());
		if (!match.hasMatch())
			return {};
		const auto date = QDate::fromString(match.captured(1), Qt::ISODate);
		const auto time = QTime::fromString(match.captured(2), "HH:mm:ss");
		if (!date.isValid() || !time.isValid())
			return {};
		const double fraction = match.captured(3).isEmpty() ? 0.0 : ("0." + match.captured(3)).toDouble();
		const double mjd = date.toJulianDay() - 2400001.0 + (time.msecsSinceStartOfDay() / 1000.0 + fraction) / 86400.0;
		const double tdb = orbital_core::time::toMjdTdb(mjd, orbital_core::time::Scale::UTC);
		return std::isfinite(tdb) ? std::optional<double>(tdb) : std::nullopt;
	}

	void writeFile(const QString& path, const QByteArray& bytes)
	{
		QSaveFile file(path);
		if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
			throw std::runtime_error(
			    (QStringLiteral("Cannot save close approaches: ") + file.errorString()).toStdString());
	}

	QJsonValue optionalNumber(std::optional<double> value)
	{
		return value && std::isfinite(*value) ? QJsonValue(*value) : QJsonValue(QJsonValue::Null);
	}

	void saveResult(const QString& directory, const QString& designation, const orbital_core::FitPipelineResult& source,
	                const hz::CloseApproachSearchRequest& request, const hz::CloseApproachSearchResult& result,
	                double seconds)
	{
		QJsonArray events;
		for (const auto& row : result.events) {
			events.append(
			    QJsonObject{{"body", QString::fromStdString(std::string(hz::longLabel(row.encounter.planet)))},
			                {"mjdTdb", row.mjdTdb},
			                {"distanceAu", row.encounter.distAu},
			                {"minimum3SigmaAu", optionalNumber(row.minimumDistanceAu)},
			                {"maximum3SigmaAu", optionalNumber(row.maximumDistanceAu)},
			                {"time3SigmaDays", optionalNumber(row.timeUncertaintyDays)},
			                {"relativeSpeedKmS", row.relativeSpeedKmS},
			                {"vInfinityKmS",
			                 optionalNumber(row.encounter.vInfKmS > 0.0 ? std::optional<double>(row.encounter.vInfKmS)
			                                                            : std::nullopt)},
			                {"uncertaintyNote", QString::fromStdString(row.uncertaintyNote)}});
		}
		const auto& fit = source.fit;
		const auto r = fit.fitEpochOrbit->getR(), v = fit.fitEpochOrbit->getV() * Constants::GAUSS_K;
		QJsonArray covariance;
		const auto& cov = fit.uncertainty.fitParameterCovariance;
		for (int i = 0; i < cov.rows(); ++i) {
			QJsonArray values;
			for (int j = 0; j < cov.cols(); ++j)
				values.append(optionalNumber(cov(i, j)));
			covariance.append(values);
		}
		QJsonObject limits;
		for (auto body : hz::kAllPlanets)
			limits[QString::fromStdString(std::string(hz::longLabel(body)))] =
			    request.maximumDistanceAu.value_or(hz::defaultCloseApproachLimitAu(body));
		const auto& forces = source.effectiveConfig.nbody;
		const QJsonObject object{
		    {"schema", "meridian-close-approaches-v2"},
		    {"object", designation},
		    {"timeScale", "TDB"},
		    {"startMjdTdb", request.startMjdTdb},
		    {"endMjdTdb", request.endMjdTdb},
		    {"reportingLimitsAu", limits},
		    {"seconds", seconds},
		    {"stoppedAtSurface", result.stoppedAtSurface},
		    {"uncertaintyMethod", "Local linear propagation of objective-scaled full fit covariance; 3 sigma"},
		    {"source", QJsonObject{{"epochMjdTdb", fit.fitEpochOrbit->epochMjd()},
		                           {"stateAuAuDay", QJsonArray{r.x, r.y, r.z, v.x, v.y, v.z}},
		                           {"fitParameterCovariance", covariance},
		                           {"covarianceObjectiveScale", fit.uncertainty.covarianceObjectiveScale},
		                           {"planetMask", forces.planetMask},
		                           {"stepDays", forces.stepSize},
		                           {"sb441", forces.useSb441},
		                           {"relativity", orbital_core::relativityModeName(forces.forceModel.relativity)},
		                           {"sunJ2", true},
		                           {"earthJ2", forces.forceModel.earthJ2},
		                           {"earthJ3", forces.forceModel.earthJ3},
		                           {"earthJ4", forces.forceModel.earthJ4},
		                           {"earthRotatingPole", forces.forceModel.rotatingEarthPole},
		                           {"targetNaifId", forces.targetNaifId},
		                           {"inactiveInteriorMassPolicy", static_cast<int>(forces.inactiveInteriorMassPolicy)},
		                           {"ngModel", static_cast<int>(fit.finalNonGrav.model())},
		                           {"ngEnabled", fit.finalNonGrav.enabled()},
		                           {"A1", fit.finalNonGrav.A1()},
		                           {"A2", fit.finalNonGrav.A2()},
		                           {"A3", fit.finalNonGrav.A3()},
		                           {"DT", fit.finalNonGrav.DT()},
		                           {"ngSolveMask", static_cast<int>(fit.finalNonGrav.fitParameterMask())}}},
		    {"events", events}};
		if (!QDir().mkpath(directory))
			throw std::runtime_error("Cannot create close-approach archive folder");
		writeFile(QDir(directory).filePath("result.json"), QJsonDocument(object).toJson());
		CloseApproachTableModel table;
		table.setRows(result.events);
		QString text =
		    designation + "\nDistance: au; velocity: km/s; date/time: TDB; uncertainties: local linear 3 sigma\n";
		text += table.toPlainText();
		writeFile(QDir(directory).filePath("encounters.txt"), text.toUtf8());
	}
}  // namespace

CloseApproachesViewModel::CloseApproachesViewModel(SolutionViewModel& solution, QObject* parent)
    : QObject(parent), m_solution(solution)
{
	connect(&solution, &SolutionViewModel::changed, this, [this] {
		++m_revision;
		if (m_worker.joinable())
			m_worker.request_stop();
		m_table.setRows({});
		m_hasResult = false;
		m_summary.clear();
		m_error.clear();
		m_archive.clear();
		emit changed();
	});
}

CloseApproachesViewModel::~CloseApproachesViewModel()
{
	if (m_worker.joinable()) {
		m_worker.request_stop();
		m_worker.join();
	}
}

bool CloseApproachesViewModel::available() const
{
	const auto& source = m_solution.snapshot();
	return source && source->fit.fitEpochOrbit && source->fit.converged && !source->fit.cancelled &&
	       source->effectiveConfig.nbody.enabled;
}

void CloseApproachesViewModel::setBlocked(bool value)
{
	if (m_blocked == value)
		return;
	m_blocked = value;
	emit changed();
}

QString CloseApproachesViewModel::validationError(const QString& startUtc, const QString& endUtc,
                                                  double maximumAu) const
{
	if (!available())
		return tr("A converged N-body solution is required.");
	const auto start = utcToTdb(startUtc), end = utcToTdb(endUtc);
	if (!start || !end)
		return tr("Enter valid UTC dates as YYYY-MM-DD HH:MM:SS.");
	if (*end <= *start)
		return tr("End UTC must be later than start UTC.");
	const double epoch = m_solution.snapshot()->fit.fitEpochOrbit->epochMjd();
	if (std::max(std::abs(*start - epoch), std::abs(*end - epoch)) > 365250.0)
		return tr("The interval must be within 1000 years of the fit epoch.");
	if (!std::isfinite(maximumAu) || maximumAu < 0)
		return tr("Enter a positive maximum distance.");
	return {};
}

void CloseApproachesViewModel::run(const QString& startUtc, const QString& endUtc, double maximumAu)
{
	if (!canRun())
		return;
	m_error = validationError(startUtc, endUtc, maximumAu);
	if (!m_error.isEmpty()) {
		emit changed();
		return;
	}
	hz::CloseApproachSearchRequest request{*utcToTdb(startUtc), *utcToTdb(endUtc),
	                                       maximumAu > 0 ? std::optional<double>(maximumAu) : std::nullopt};
	auto source = *m_solution.snapshot();
	source.fit.fitEpochOrbit = std::make_shared<Orbit>(*source.fit.fitEpochOrbit);
	if (source.fit.orbit)
		source.fit.orbit = std::make_shared<Orbit>(*source.fit.orbit);
	const auto revision = m_revision;
	const auto designation = m_solution.objectTitle();
	const auto folder =
	    QDir(m_solution.directory().isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)
	                                          : m_solution.directory())
	        .filePath("CloseApproaches/run-" + QUuid::createUuid().toString(QUuid::WithoutBraces));
	m_busy = true;
	m_cancelling = false;
	m_progress = tr("Preparing close-approach search...");
	emit changed();
	m_worker = std::jthread([this, source = std::move(source), request, revision, designation, folder,
	                         interval = startUtc + " to " + endUtc + " UTC"](std::stop_token stop) {
		const auto started = std::chrono::steady_clock::now();
		auto result = std::make_shared<hz::CloseApproachSearchResult>();
		QString error, archive;
		bool succeeded = false;
		try {
			orbital_core::FitRunControl control{stop, [this, revision](const std::string& status) {
				                                    QMetaObject::invokeMethod(
				                                        this,
				                                        [this, revision, status = QString::fromStdString(status)] {
					                                        if (revision != m_revision || !m_busy || m_cancelling)
						                                        return;
					                                        m_progress = status;
					                                        emit changed();
				                                        },
				                                        Qt::QueuedConnection);
			                                    }};
			*result = hz::searchCloseApproaches(source, request, control);
			if (!result->cancelled && !stop.stop_requested()) {
				succeeded = true;
				saveResult(folder, designation, source, request, *result,
				           std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
				archive = folder;
			}
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
		}
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
		QMetaObject::invokeMethod(
		    this,
		    [this, revision, result, error, archive, seconds, interval, succeeded,
		     cancelled = stop.stop_requested() || result->cancelled] {
			    m_busy = false;
			    m_cancelling = false;
			    if (revision == m_revision) {
				    m_error = cancelled ? QString() : error;
				    m_progress = cancelled   ? tr("Search cancelled")
				                 : succeeded ? tr("Search completed")
				                             : tr("Search failed");
				    if (succeeded && !cancelled) {
					    m_table.setRows(result->events);
					    m_hasResult = true;
					    m_archive = archive;
					    m_summary = tr("%1 encounters | %2 s | %3")
					                    .arg(result->events.size())
					                    .arg(seconds, 0, 'f', 1)
					                    .arg(interval);
					    if (result->stoppedAtSurface)
						    m_summary += tr(" | Nominal surface crossing; later encounters omitted");
					    emit completed();
				    }
			    }
			    emit changed();
		    },
		    Qt::QueuedConnection);
	});
}

void CloseApproachesViewModel::cancel()
{
	if (!m_busy || m_cancelling)
		return;
	m_cancelling = true;
	m_progress = tr("Cancelling...");
	m_worker.request_stop();
	emit changed();
}

void CloseApproachesViewModel::openArchive()
{
	if (!m_archive.isEmpty() && !QDesktopServices::openUrl(QUrl::fromLocalFile(m_archive))) {
		m_error = tr("Cannot open the close-approach archive.");
		emit changed();
	}
}
