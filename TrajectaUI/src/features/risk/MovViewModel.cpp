#include "MovViewModel.h"
#include <orbital_core/hazard/ManifoldOfVariations.h>
#include <orbital_core/TimeScale.h>
#include <QCoreApplication>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QFileInfo>
#include <QStandardPaths>
#include <QUuid>
#include <QRegularExpression>
#include <QTimeZone>
#include <chrono>
#include <cmath>

namespace {
	namespace oc = orbital_core;
	namespace hz = oc::hazard;

	QString utc(double tdb)
	{
		const double days = oc::time::fromMjdTdb(tdb, oc::time::Scale::UTC) - 40587;
		return QDateTime::fromSecsSinceEpoch(static_cast<qint64>(std::ceil(days * 86400)), QTimeZone::UTC)
		    .toString("yyyy-MM-dd HH:mm:ss");
	}

	double epoch(const QString& text)
	{
		static const QRegularExpression pattern(
		    "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}(?:\\.[0-9]{1,3})?$");
		if (!pattern.match(text.trimmed()).hasMatch())
			throw std::invalid_argument("Enter UTC as YYYY-MM-DD HH:MM:SS");
		const auto date = QDateTime::fromString(text.trimmed().replace(' ', 'T') + "Z", Qt::ISODateWithMs);
		if (!date.isValid())
			throw std::invalid_argument("Invalid UTC date");
		const double mjd = date.date().toJulianDay() - 2400001.0 + date.time().msecsSinceStartOfDay() / 86400000.0;
		return oc::time::toMjdTdb(mjd, oc::time::Scale::UTC);
	}

	int integer(const QVariantMap& options, const char* key, int min, int max)
	{
		bool ok = false;
		const auto d = options.value(key).toDouble(&ok);
		if (!ok || !std::isfinite(d) || d != std::floor(d) || d < min || d > max)
			throw std::invalid_argument(std::string("Invalid MOV option: ") + key);
		return static_cast<int>(d);
	}
}  // namespace

MovViewModel::MovViewModel(SolutionViewModel& solution, std::function<OptimizationSettings()> settings, QObject* parent)
    : QObject(parent), m_solution(solution), m_settings(std::move(settings))
{
	connect(&solution, &SolutionViewModel::changed, this, [this] {
		if (!m_external)
			resetSource();
	});
}

MovViewModel::~MovViewModel()
{
	if (m_worker.joinable()) {
		m_worker.request_stop();
		m_worker.join();
	}
}

void MovViewModel::resetSource()
{
	++m_revision;
	if (m_worker.joinable())
		m_worker.request_stop();
	m_runs.clear();
	m_selected = -1;
	m_error.clear();
	emit changed();
}

QString MovViewModel::sourceName() const
{
	return m_external ? QFileInfo(m_file).fileName() : m_solution.objectTitle();
}

int MovViewModel::observationCount() const
{
	return m_external ? static_cast<int>(m_external->all.size()) : m_solution.total();
}

QString MovViewModel::seedUtc() const
{
	const auto* obs = m_external              ? &m_external->all
	                  : m_solution.snapshot() ? &m_solution.snapshot()->observations.all
	                                          : nullptr;
	if (!obs || obs->empty())
		return {};
	try {
		return utc(obs->back().epochMjdTdb());
	} catch (...) {
		return {};
	}
}

QVariantMap MovViewModel::result() const
{
	return m_selected >= 0 && m_selected < m_runs.size() ? m_runs[m_selected].toMap() : QVariantMap{};
}

void MovViewModel::setSelectedRun(int i)
{
	if (i >= 0 && i < m_runs.size() && i != m_selected) {
		m_selected = i;
		emit changed();
	}
}

void MovViewModel::setBlocked(bool b)
{
	if (m_blocked != b) {
		m_blocked = b;
		emit changed();
	}
}

void MovViewModel::selectFile(const QUrl& url)
{
	if (m_busy || m_blocked)
		return;
	const QFileInfo file(url.toLocalFile());
	if (!url.isLocalFile() || !file.isReadable() || !file.isFile()) {
		m_error = tr("Cannot read observations");
		emit changed();
		return;
	}
	m_busy = true;
	m_cancelling = false;
	m_error.clear();
	m_progress = tr("Loading MOV observations...");
	emit changed();
	const auto path = file.absoluteFilePath();
	m_worker = std::jthread([this, path](std::stop_token stop) {
		auto observations = std::make_shared<oc::ObservationSet>();
		QString error;
		try {
			*observations = oc::loadObservations(std::filesystem::path(path.toStdWString()));
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
		}
		QMetaObject::invokeMethod(
		    this,
		    [this, observations, path, error, cancelled = stop.stop_requested()] {
			    m_busy = false;
			    m_cancelling = false;
			    if (!cancelled && error.isEmpty()) {
				    m_external = std::move(*observations);
				    m_file = path;
				    resetSource();
				    m_progress.clear();
			    } else {
				    m_error = cancelled ? QString() : error;
				    m_progress = cancelled ? tr("Cancelled") : tr("Loading failed");
			    }
			    emit changed();
		    },
		    Qt::QueuedConnection);
	});
}

void MovViewModel::useAcceptedSolution()
{
	if (m_busy || m_blocked || !m_external)
		return;
	m_external.reset();
	m_file.clear();
	resetSource();
}

void MovViewModel::run(const QVariantMap& options)
{
	if (!canRun())
		return;
	hz::MovConfig config;
	auto source = std::make_shared<oc::ObservationSet>();
	oc::FitConfig settings;
	std::optional<oc::FitResult> nominal;
	bool prepare = false, refitNominal = false;
	OptimizationSettings rawSettings;
	try {
		config.startMjdTdb = epoch(options.value("startUtc").toString());
		config.endMjdTdb = epoch(options.value("endUtc").toString());
		if (config.endMjdTdb <= config.startMjdTdb || config.endMjdTdb - config.startMjdTdb > 30)
			throw std::invalid_argument("MOV interval must be positive and at most 30 days");
		config.rangeNodes = integer(options, "rangeNodes", 8, 250);
		config.rateNodes = integer(options, "rateNodes", 8, 250);
		config.workers = integer(options, "workers", 1, 64);
		config.sampling = static_cast<hz::MovSampling>(integer(options, "sampling", 0, 2));
		const int target = integer(options, "target", 0, 2);
		config.target = target == 0 ? hz::Planet::Earth : target == 1 ? hz::Planet::Moon : hz::Planet::Mars;
		config.altitudeKm = target == 0 && integer(options, "boundary", 0, 1) == 1 ? 100 : 0;
		const bool all = integer(options, "observationSelection", 0, 1) == 1;
		if (m_external) {
			*source = *m_external;
			prepare = true;
			rawSettings = m_settings();
			settings = OrbitWorker::makeFitConfig(rawSettings);
			settings.nbody.targetNaifId = source->targetNaifId;
		} else {
			const auto& accepted = *m_solution.snapshot();
			settings = accepted.effectiveConfig;
			*source = accepted.observations;
			if (!all) {
				source->all.clear();
				for (const auto& o : accepted.fittedObservations())
					source->all.push_back(o);
			}
			if (config.sampling != hz::MovSampling::Grid) {
				nominal = accepted.fit;
				if (nominal->fitEpochOrbit)
					nominal->fitEpochOrbit = std::make_shared<Orbit>(*nominal->fitEpochOrbit);
				refitNominal = all && accepted.fittedObservations().size() != source->all.size();
			}
		}
		if (settings.nonGrav.enabled || !settings.nbody.enabled || !settings.useWeights)
			throw std::invalid_argument(
			    "MOV requires N-body dynamics, weighted observations and non-gravitational parameters disabled");
		if (source->all.size() < 3)
			throw std::invalid_argument("MOV needs at least three selected optical observations");
	} catch (const std::exception& e) {
		m_error = QString::fromUtf8(e.what());
		emit changed();
		return;
	}
	const auto revision = m_revision;
	const auto name = sourceName();
	const auto root = m_external || m_solution.directory().isEmpty()
	                      ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation)
	                      : m_solution.directory();
	const auto directory = QDir(root).filePath("MOV/run-" + QUuid::createUuid().toString(QUuid::WithoutBraces));
	const auto runtime = QCoreApplication::applicationDirPath().toStdString();
	m_busy = true;
	m_cancelling = false;
	m_error.clear();
	m_progress = refitNominal ? tr("Fitting selected MOV observations...") : tr("Constructing admissible region...");
	emit changed();
	emit runStarted();
	m_worker = std::jthread([this, source, settings, nominal = std::move(nominal), config, prepare, refitNominal,
	                         rawSettings, revision, name, directory, runtime](std::stop_token stop) mutable {
		const auto start = std::chrono::steady_clock::now();
		QVariantMap summary;
		QString error;
		auto lastUpdate = start;
		try {
			if (prepare)
				*source = oc::prepareObservationCopy(*source, {.debias = rawSettings.debiasEnabled,
				                                               .weightProfile = settings.weightProfile,
				                                               .runtimeDirectory = runtime});
			if (refitNominal && nominal->fitEpochOrbit && !stop.stop_requested()) {
				std::vector<oc::ObservationOverride> included;
				included.reserve(source->all.size());
				for (const auto& o : source->all)
					included.push_back({o.id, oc::ObservationDisposition::Include});
				oc::FitRequest request{nominal->fitEpochOrbit, source->all};
				request.options.selection = oc::FitSelectionMode::Fixed;
				request.observationOverrides = included;
				nominal = oc::FitSession{settings}.optimize(request, {.cancel = stop});
			}
			const auto r = hz::computeManifoldOfVariations(
			    source->all, settings, config, nominal ? &*nominal : nullptr, stop,
			    [this, revision, &lastUpdate](int done, int total) {
				    const auto now = std::chrono::steady_clock::now();
				    if (done < total && now - lastUpdate < std::chrono::milliseconds(150))
					    return;
				    lastUpdate = now;
				    QMetaObject::invokeMethod(
				        this,
				        [this, revision, done, total] {
					        if (m_revision == revision && m_busy && !m_cancelling) {
						        m_progress = tr("MOV nodes: %1 / %2").arg(done).arg(total);
						        emit changed();
					        }
				        },
				        Qt::QueuedConnection);
			    });
			int fits = 0, selected = 0, hits = 0, failures = 0;
			for (const auto& n : r.nodes) {
				fits += std::isfinite(n.chiSquare);
				selected += n.weight > 0;
				hits += n.status == hz::MovNodeStatus::Complete && n.entry.has_value();
				failures +=
				    n.status == hz::MovNodeStatus::FitFailed || n.status == hz::MovNodeStatus::PropagationFailed;
			}
			QStringList warnings;
			for (const auto& w : r.warnings)
				warnings.append(QString::fromStdString(w));
			summary = {
			    {"label", QDateTime::currentDateTimeUtc().toString("yyyy-MM-dd HH:mm:ss 'UTC'") + " | " + name},
			    {"sampling", r.sampling == hz::MovSampling::Cobweb ? "Cobweb" : "Grid"},
			    {"nodes", static_cast<int>(r.nodes.size())},
			    {"fits", fits},
			    {"selected", selected},
			    {"hits", hits},
			    {"failures", failures},
			    {"cancelled", r.cancelled},
			    {"warnings", warnings.join('\n')},
			    {"curvature", r.curvatureSnr},
			    {"probability", r.probabilityLower ? QString::number(*r.probabilityLower, 'g', 9) : tr("Unavailable")},
			    {"probabilityUpper",
			     r.probabilityUpper ? QString::number(*r.probabilityUpper, 'g', 9) : tr("Unavailable")},
			    {"seconds", std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()}};
			hz::saveManifoldOfVariations(r, source->all, settings, std::filesystem::path(directory.toStdWString()),
			                             name.toStdString());
			summary["archive"] = directory;
			hz::renderManifoldOfVariationsReport(std::filesystem::path(directory.toStdWString()));
			summary["report"] = QDir(directory).filePath("report.html");
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
		}
		QMetaObject::invokeMethod(
		    this,
		    [this, revision, summary, error, cancelled = stop.stop_requested()] {
			    m_busy = false;
			    m_cancelling = false;
			    if (revision == m_revision) {
				    m_error = error;
				    m_progress = cancelled ? tr("Cancelled") : error.isEmpty() ? tr("MOV completed") : tr("MOV error");
				    if (!summary.isEmpty()) {
					    m_runs.append(summary);
					    m_selected = m_runs.size() - 1;
				    }
			    }
			    emit changed();
			    if (revision == m_revision && summary.isEmpty() && !error.isEmpty())
				    emit runFailed();
		    },
		    Qt::QueuedConnection);
	});
}

void MovViewModel::cancel()
{
	if (m_busy && !m_cancelling) {
		m_cancelling = true;
		m_worker.request_stop();
		m_progress = tr("Cancelling...");
		emit changed();
	}
}

void MovViewModel::openReport()
{
	const auto path = result().value("report").toString();
	if (!path.isEmpty() && !QDesktopServices::openUrl(QUrl::fromLocalFile(path))) {
		m_error = tr("Cannot open MOV report");
		emit changed();
	}
}

void MovViewModel::openArchive()
{
	const auto path = result().value("archive").toString();
	if (!path.isEmpty() && !QDesktopServices::openUrl(QUrl::fromLocalFile(path))) {
		m_error = tr("Cannot open MOV archive");
		emit changed();
	}
}
