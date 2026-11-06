#include "AtmosphereViewModel.h"

#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QStandardPaths>
#include <QUuid>
#include <cmath>
#include <numbers>

AtmosphereViewModel::~AtmosphereViewModel()
{
	m_worker.request_stop();
	if (m_worker.joinable())
		m_worker.join();
}

void AtmosphereViewModel::setBlocked(bool value)
{
	if (value == m_blocked)
		return;
	m_blocked = value;
	emit changed();
}

QStringList AtmosphereViewModel::events() const
{
	QStringList labels;
	if (m_archive)
		for (const auto& e : m_archive->events)
			labels.push_back(tr("%1 to %2 UTC | %3 entries")
			                     .arg(atmosphere::formatTdb(e.startMjdTdb), atmosphere::formatTdb(e.endMjdTdb))
			                     .arg(e.entries.size()));
	return labels;
}

void AtmosphereViewModel::setSelectedEvent(int value)
{
	if (busy() || !m_archive || value < 0 || value >= m_archive->events.size() || value == m_selectedEvent)
		return;
	m_selectedEvent = value;
	emit changed();
}

QString AtmosphereViewModel::sourceSummary() const
{
	if (m_sourceMode != 0 || !m_archive || m_selectedEvent < 0)
		return {};
	return tr("%1 entry states | ECEF / UTC | H: %2")
	    .arg(m_archive->events[m_selectedEvent].entries.size())
	    .arg(m_magnitudeEdited      ? tr("manual")
	         : m_archive->magnitude ? tr("archived solution")
	                                : tr("not available"));
}

QStringList AtmosphereViewModel::runs() const
{
	QStringList labels;
	for (const auto& run : m_runs)
		labels.push_back(run.label);
	return labels;
}

void AtmosphereViewModel::setSelectedRun(int value)
{
	if (value < 0 || value >= m_runs.size() || value == m_selectedRun)
		return;
	m_selectedRun = value;
	emit changed();
}

QString AtmosphereViewModel::summary() const { return m_selectedRun >= 0 ? m_runs[m_selectedRun].summary : QString(); }

QString AtmosphereViewModel::runError() const { return m_selectedRun >= 0 ? m_runs[m_selectedRun].error : QString(); }

QString AtmosphereViewModel::resultPath() const { return m_selectedRun >= 0 ? m_runs[m_selectedRun].path : QString(); }

bool AtmosphereViewModel::reportAvailable() const
{
	return QFileInfo(resultPath() + "/results/maps/field.html").isFile();
}

void AtmosphereViewModel::openResult(bool report)
{
	if (resultPath().isEmpty() || (report && !reportAvailable()))
		return;
	const auto path = report ? resultPath() + "/results/maps/field.html" : resultPath();
	if (!QDesktopServices::openUrl(QUrl::fromLocalFile(path))) {
		m_error = tr("Cannot open atmospheric results.");
		emit changed();
	}
}

void AtmosphereViewModel::cancel()
{
	if (!busy() || m_cancelling)
		return;
	m_cancelling = true;
	m_progress = tr("Cancelling...");
	m_worker.request_stop();
	emit changed();
}

void AtmosphereViewModel::run()
{
	if (!canRun())
		return;
	strewn::JobSpec job;
	std::vector<strewn::EntryState> entries;
	QJsonObject provenance;
	const auto source = m_sourceMode == 0 ? m_archive : nullptr;
	QString root;
	try {
		const auto parameters = parameterOptions();
		job.physics = parameters.physics;
		job.integration = parameters.integration;
		job.monteCarlo = parameters.monteCarlo;
		job.maps = parameters.maps;
		std::visit([&job](const auto& environment) { job.environment = environment; }, parameters.environment);
		if (source) {
			const auto& event = source->events.at(m_selectedEvent);
			entries = event.entries;
			provenance = source->provenance;
			provenance["selected_event"] = m_selectedEvent;
			provenance["event_start_mjd_tdb"] = event.startMjdTdb;
			provenance["event_end_mjd_tdb"] = event.endMjdTdb;
			provenance["selected_entries"] = qint64(entries.size());
			auto directory = QFileInfo(source->path).absoluteDir();
			auto parent = directory;
			if (parent.cdUp() && parent.dirName() == "MC" && parent.cdUp())
				directory = parent;
			root = directory.filePath("Atmosphere");
			job.maps->title = QFileInfo(source->path).absoluteDir().dirName().toStdString();
		} else {
			const auto n = [this](const char* key) { return m_values.value(key).toString().toDouble(); };
			constexpr double rad = std::numbers::pi / 180;
			const strewn::GeodeticPosition geo{n("latitude") * rad, n("longitude") * rad, 100000};
			const double speed = n("speed") * 1000, heading = n("heading") * rad, slope = n("slope") * rad;
			const auto velocity =
			    strewn::enuToEcef({speed * std::cos(slope) * std::sin(heading),
			                       speed * std::cos(slope) * std::cos(heading), speed * std::sin(slope)},
			                      geo.latitudeRad, geo.longitudeRad);
			entries.push_back(
			    {strewn::ecefPosition(geo), velocity, atmosphere::parseUtc(m_values["epoch"].toString()), 1, 0});
			atmosphere::validateEntry(entries.front());
			provenance = {{"schema", "atmospheric-fixed-entry-1"},
			              {"conditioning", "user supplied fixed entry"},
			              {"frame", "ECEF"},
			              {"epoch_scale", "MJD UTC"},
			              {"height_wgs84_m", 100000}};
			root = QDir(m_solution.directory().isEmpty()
			                ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)
			                : m_solution.directory())
			           .filePath("Atmosphere");
			job.maps->title = "Fixed atmospheric entry";
		}
		provenance["size_source"] = m_sizeMode != 0                ? "user supplied size range"
		                            : source && !m_magnitudeEdited ? "archived source OEF H"
		                                                           : "user supplied H";
		provenance["values"] = QJsonObject::fromVariantMap(m_values);
		job.environmentDirectory =
		    std::filesystem::path(QDir(QStandardPaths::writableLocation(QStandardPaths::CacheLocation))
		                              .filePath("StrewnField")
		                              .toStdWString());
	} catch (const std::exception& e) {
		m_error = QString::fromUtf8(e.what());
		emit changed();
		return;
	}
	const auto label = QDateTime::currentDateTimeUtc().toString("yyyy-MM-dd HH:mm:ss 'UTC'") +
	                   tr(" | Atmospheric MC | N=%1").arg(job.monteCarlo.runs);
	m_busy = true;
	m_cancelling = false;
	m_progressFraction = 0;
	m_progress = tr("Preparing atmospheric entries and environment...");
	m_error.clear();
	emit changed();
	emit runStarted();
	const auto revision = m_sourceRevision;
	m_worker = std::jthread([this, job = std::move(job), entries = std::move(entries), provenance, root, label, source,
	                         revision](std::stop_token stop) mutable {
		QString path, error, summary;
		try {
			if (source)
				atmosphere::verifyArchive(*source);
			if (stop.stop_requested())
				throw std::runtime_error("Cancelled before atmospheric preparation");
			if (job.monteCarlo.convergence)
				job.monteCarlo.convergence = strewn::atmosphericConvergenceGrid(entries, *job.monteCarlo.convergence);
			path = QDir(root).filePath("run-" + QUuid::createUuid().toString(QUuid::WithoutBraces));
			if (!QDir().mkpath(path))
				throw std::runtime_error("Cannot create atmospheric output directory");
			atmosphere::saveEntries(path + "/entries.json", entries);
			atmosphere::saveJson(path + "/source.json", provenance);
			job.input = strewn::EntryFile{std::filesystem::path((path + "/entries.json").toStdWString())};
			job.outputDirectory = std::filesystem::path((path + "/results").toStdWString());
			if (auto* environment = std::get_if<strewn::EnvironmentFile>(&job.environment)) {
				const auto target = path + "/environment-input.json";
				if (!QFile::copy(QString::fromStdWString(environment->path.wstring()), target))
					throw std::runtime_error("Cannot snapshot the selected environment");
				environment->path = std::filesystem::path(target.toStdWString());
			}
			strewn::saveJobConfiguration(std::filesystem::path((path + "/launch.json").toStdWString()), job);
			strewn::RunControl control;
			control.cancel = stop;
			auto last = std::chrono::steady_clock::now();
			control.progress = [this, &last, revision](const strewn::RunSummary& s) {
				const auto now = std::chrono::steady_clock::now();
				if (s.completed != s.requested && !s.converged && now - last < std::chrono::milliseconds(100))
					return;
				last = now;
				const auto text = tr("Atmospheric MC: %1 / %2 | %3 ground impacts | %4 failed | %5 incomplete")
				                      .arg(s.completed)
				                      .arg(s.requested)
				                      .arg(s.trackedImpacts)
				                      .arg(s.failed)
				                      .arg(s.truncated);
				const double fraction = s.requested ? double(s.completed) / s.requested : 0;
				QMetaObject::invokeMethod(
				    this,
				    [this, text, fraction, revision] {
					    if (revision != m_sourceRevision)
						    return;
					    if (!m_cancelling)
						    m_progress = text;
					    m_progressFraction = fraction;
					    emit changed();
				    },
				    Qt::QueuedConnection);
			};
			const auto result = strewn::StrewnJob{}.run(job, control);
			const auto& s = result.simulation;
			summary = tr("%1 / %2 realizations | %3 with ground impacts | %4 tracked impacts | %5 failed | %6 "
			             "incomplete | %7 s")
			              .arg(s.completed)
			              .arg(s.requested)
			              .arg(s.runsWithImpacts)
			              .arg(s.trackedImpacts)
			              .arg(s.failed)
			              .arg(s.truncated)
			              .arg(s.elapsedS, 0, 'f', 1);
			summary += s.cancelled   ? tr(" | Cancelled")
			           : s.converged ? tr(" | JS converged")
			                         : tr(" | Run limit reached");
			if (s.convergenceBlocked)
				summary += tr(" | JS blocked by incomplete/out-of-grid results");
			if (s.wallTimeouts)
				summary += tr(" | %1 timeouts").arg(s.wallTimeouts);
			error = QString::fromStdString(result.artifacts.error);
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
			summary = stop.stop_requested() ? tr("Cancelled") : tr("Failed");
		}
		if (!path.isEmpty())
			try {
				atmosphere::saveJson(path + "/outcome.json",
				                     {{"label", label}, {"summary", summary}, {"error", error}});
			} catch (const std::exception& e) {
				error += "\n" + QString::fromUtf8(e.what());
			}
		QMetaObject::invokeMethod(
		    this,
		    [this, path, label, summary, error, revision] {
			    m_busy = m_cancelling = false;
			    if (revision != m_sourceRevision) {
				    emit changed();
				    return;
			    }
			    m_error = error;
			    m_progress = summary;
			    m_runs.push_back({label, path, summary, error});
			    m_selectedRun = static_cast<int>(m_runs.size()) - 1;
			    emit changed();
		    },
		    Qt::QueuedConnection);
	});
}
