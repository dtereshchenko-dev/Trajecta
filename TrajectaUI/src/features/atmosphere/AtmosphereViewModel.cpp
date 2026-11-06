#include "AtmosphereViewModel.h"

#include <QDateTime>
#include <QFileInfo>
#include <QRegularExpression>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
	struct RangeField {
		const char* key;
		strewn::Range strewn::ParameterRanges::* member;
		double scale;
	};

	constexpr RangeField ranges[]{{"density", &strewn::ParameterRanges::densityKgM3, 1},
	                              {"strength", &strewn::ParameterRanges::strongStrengthPa, 1e6},
	                              {"ablation", &strewn::ParameterRanges::ablationKgJ, 1e-8},
	                              {"dust", &strewn::ParameterRanges::dustFraction, 1},
	                              {"largest", &strewn::ParameterRanges::largestFraction, 1},
	                              {"beta", &strewn::ParameterRanges::beta, 1},
	                              {"wind", &strewn::ParameterRanges::windScale, 1},
	                              {"drag", &strewn::ParameterRanges::dragScale, 1}};

	double number(const QVariantMap& values, const QString& key)
	{
		bool ok;
		const auto value = values.value(key).toString().toDouble(&ok);
		if (!ok || !std::isfinite(value))
			throw std::invalid_argument(("Invalid number: " + key).toStdString());
		return value;
	}

	std::uint64_t integer(const QVariantMap& values, const QString& key, std::uint64_t maximum)
	{
		const auto text = values.value(key).toString();
		bool ok;
		const auto value = text.toULongLong(&ok);
		if (!ok || !QRegularExpression(QStringLiteral("^[0-9]+$")).match(text).hasMatch() || value > maximum)
			throw std::invalid_argument(("Invalid integer: " + key).toStdString());
		return value;
	}

	QString localPath(const QUrl& url)
	{
		return url.isLocalFile() ? QFileInfo(url.toLocalFile()).absoluteFilePath() : QString();
	}
}  // namespace

AtmosphereViewModel::AtmosphereViewModel(SolutionViewModel& solution, RiskMonitorViewModel& risk, QObject* parent)
    : QObject(parent), m_solution(solution), m_risk(risk)
{
	m_values = {{"epoch", ""},      {"latitude", "0"},        {"longitude", "0"},    {"speed", "12"},
	            {"heading", "90"},  {"slope", "-30"},         {"magnitude", ""},     {"massMin", "100"},
	            {"massMax", "100"}, {"diameterMin", "1"},     {"diameterMax", "1"},  {"alpha", "0.2"},
	            {"minMass", "1"},   {"representatives", "3"}, {"runs", "1000"},      {"workers", "15"},
	            {"timeout", "300"}, {"seed", "42"},           {"threshold", "0.02"}, {"batch", "50"},
	            {"checks", "3"},    {"cell", "25"},           {"sigma", "60"},       {"groundHeight", "0"}};
	strewn::JobSpec job;
	job.monteCarlo.parameters.size = strewn::PhotometricSize{};
	strewn::applyMaterialPreset(job, strewn::MaterialPreset::OrdinaryChondrite);
	readMaterial(job);
	connect(&m_solution, &SolutionViewModel::changed, this, &AtmosphereViewModel::refreshSource);
	connect(&m_risk, &RiskMonitorViewModel::changed, this, [this] {
		const auto sources = m_risk.atmosphericSources();
		if (sources == m_sources)
			return;
		m_sources = sources;
		emit changed();
	});
	m_sources = m_risk.atmosphericSources();
	refreshSource();
}

void AtmosphereViewModel::refreshSource()
{
	const auto designation = m_solution.designation();
	if (!designation.isEmpty()) {
		if (!m_sourceObject.isEmpty() && designation != m_sourceObject) {
			++m_sourceRevision;
			if (m_worker.joinable())
				m_worker.request_stop();
			m_archive.reset();
			m_archivePath.clear();
			m_archiveError.clear();
			m_selectedEvent = -1;
			m_runs.clear();
			m_selectedRun = -1;
			m_error.clear();
			m_progress.clear();
			emit changed();
		}
		m_sourceObject = designation;
	}
	refreshMagnitude();
}

void AtmosphereViewModel::refreshMagnitude()
{
	if (m_magnitudeEdited || busy() || m_sourceMode == 0)
		return;
	const auto magnitude = m_solution.absoluteMagnitude();
	const QString value =
	    m_solution.available() && std::isfinite(magnitude) ? QString::number(magnitude, 'g', 12) : QString();
	if (m_values.value("magnitude").toString() == value)
		return;
	m_values["magnitude"] = value;
	emit changed();
}

void AtmosphereViewModel::setValue(const QString& key, const QString& value)
{
	if (busy() || !m_values.contains(key) || m_values[key].toString() == value)
		return;
	m_values[key] = value;
	if (key == "magnitude")
		m_magnitudeEdited = true;
	emit changed();
}

void AtmosphereViewModel::readMaterial(const strewn::JobSpec& job)
{
	for (const auto& field : ranges) {
		const auto range = job.monteCarlo.parameters.*field.member;
		m_values[QString::fromLatin1(field.key) + "Min"] = QString::number(range.min / field.scale, 'g', 12);
		m_values[QString::fromLatin1(field.key) + "Max"] = QString::number(range.max / field.scale, 'g', 12);
	}
	const auto& size = std::get<strewn::PhotometricSize>(job.monteCarlo.parameters.size);
	m_values["albedoMin"] = QString::number(size.albedoMin, 'g', 12);
	m_values["albedoMode"] = QString::number(size.albedoMode, 'g', 12);
	m_values["albedoMax"] = QString::number(size.albedoMax, 'g', 12);
}

void AtmosphereViewModel::selectPreset(int index)
{
	if (busy() || index < 0 || index > 2)
		return;
	strewn::JobSpec job;
	job.monteCarlo.parameters.size = strewn::PhotometricSize{};
	strewn::applyMaterialPreset(job, static_cast<strewn::MaterialPreset>(index));
	readMaterial(job);
	emit changed();
}

int AtmosphereViewModel::preset() const
{
	try {
		strewn::JobSpec job;
		for (const auto& f : ranges)
			job.monteCarlo.parameters.*f.member = {number(m_values, QString(f.key) + "Min") * f.scale,
			                                       number(m_values, QString(f.key) + "Max") * f.scale};
		job.monteCarlo.parameters.size = strewn::PhotometricSize{
		    0, number(m_values, "albedoMin"), number(m_values, "albedoMode"), number(m_values, "albedoMax")};
		return static_cast<int>(strewn::materialPreset(job));
	} catch (...) {
		return 3;
	}
}

void AtmosphereViewModel::setSourceMode(int value)
{
	if (busy() || value < 0 || value > 1 || value == m_sourceMode)
		return;
	m_sourceMode = value;
	if (value == 0 && !m_magnitudeEdited)
		m_values["magnitude"] =
		    m_archive && m_archive->magnitude ? QString::number(*m_archive->magnitude, 'g', 15) : QString();
	refreshMagnitude();
	emit changed();
}

void AtmosphereViewModel::setSizeMode(int value)
{
	if (busy() || value < 0 || value > 2 || value == m_sizeMode)
		return;
	m_sizeMode = value;
	emit changed();
}

void AtmosphereViewModel::setEnvironmentMode(int value)
{
	if (busy() || value < 0 || value > 2 || value == m_environmentMode)
		return;
	m_environmentMode = value;
	emit changed();
}

void AtmosphereViewModel::setOffline(bool value)
{
	if (busy() || value == m_offline)
		return;
	m_offline = value;
	emit changed();
}

void AtmosphereViewModel::setConvergence(bool value)
{
	if (busy() || value == m_convergence)
		return;
	m_convergence = value;
	emit changed();
}

void AtmosphereViewModel::setMapMethod(int value)
{
	if (busy() || value < 0 || value > 1 || value == m_mapMethod)
		return;
	m_mapMethod = value;
	emit changed();
}

void AtmosphereViewModel::selectArchive(const QUrl& url)
{
	const auto path = localPath(url);
	if (busy() || path.isEmpty())
		return;
	m_archivePath = path;
	m_archive.reset();
	m_archiveError.clear();
	m_error.clear();
	m_selectedEvent = -1;
	m_values["magnitude"] = QString();
	m_magnitudeEdited = false;
	m_loading = true;
	m_cancelling = false;
	m_progressFraction = 0;
	m_progress = tr("Reading orbital MC entries...");
	emit changed();
	const auto revision = m_sourceRevision;
	m_worker = std::jthread([this, path, revision](std::stop_token stop) {
		std::shared_ptr<const atmosphere::EntryArchive> archive;
		QString error;
		try {
			archive = std::make_shared<atmosphere::EntryArchive>(atmosphere::readArchive(path, stop));
		} catch (const std::exception& e) {
			error = QString::fromUtf8(e.what());
		}
		const bool cancelled = stop.stop_requested();
		QMetaObject::invokeMethod(
		    this,
		    [this, archive, error, cancelled, revision] {
			    m_loading = m_cancelling = false;
			    if (revision != m_sourceRevision) {
				    emit changed();
				    return;
			    }
			    m_archive = cancelled ? nullptr : archive;
			    m_archiveError = cancelled ? tr("Archive loading cancelled. Select the archive again.") : error;
			    m_progress = m_archive ? tr("Entry archive ready") : m_archiveError;
			    if (m_archive) {
				    m_selectedEvent = m_archive->events.size() == 1 ? 0 : -1;
				    if (m_archive->magnitude)
					    m_values["magnitude"] = QString::number(*m_archive->magnitude, 'g', 15);
			    }
			    emit changed();
		    },
		    Qt::QueuedConnection);
	});
}

void AtmosphereViewModel::selectEnvironment(const QUrl& url)
{
	const auto path = localPath(url);
	if (busy() || path.isEmpty() || path == m_environmentPath)
		return;
	m_environmentPath = path;
	emit changed();
}

AtmosphereParameters AtmosphereViewModel::parameterOptions() const
{
	AtmosphereParameters job;
	auto& mc = job.monteCarlo;
	for (const auto& f : ranges)
		mc.parameters.*f.member = {number(m_values, QString(f.key) + "Min") * f.scale,
		                           number(m_values, QString(f.key) + "Max") * f.scale};
	if (m_sizeMode == 0)
		mc.parameters.size = strewn::PhotometricSize{number(m_values, "magnitude"), number(m_values, "albedoMin"),
		                                             number(m_values, "albedoMode"), number(m_values, "albedoMax")};
	else if (m_sizeMode == 1)
		mc.parameters.size = strewn::MassSize{{number(m_values, "massMin"), number(m_values, "massMax")}};
	else
		mc.parameters.size = strewn::DiameterSize{{number(m_values, "diameterMin"), number(m_values, "diameterMax")}};
	mc.parameters.weakProbability = 0;
	mc.runs = integer(m_values, "runs", 1000000);
	mc.workers = static_cast<unsigned>(integer(m_values, "workers", 1024));
	if (mc.runs == 0 || mc.workers == 0)
		throw std::invalid_argument("Realizations and workers must be positive");
	mc.seed = integer(m_values, "seed", std::numeric_limits<std::uint64_t>::max());
	mc.runTimeoutS = number(m_values, "timeout");
	job.physics.weibullAlpha = number(m_values, "alpha");
	job.physics.minMassKg = number(m_values, "minMass") / 1000;
	job.physics.representativesPerBin = integer(m_values, "representatives", 1000);
	if (m_convergence) {
		strewn::ConvergenceConfig convergence;
		convergence.threshold = number(m_values, "threshold");
		convergence.batchSize = integer(m_values, "batch", 1000000);
		convergence.consecutiveChecks = integer(m_values, "checks", 1000000);
		mc.convergence = convergence;
	}
	job.maps.method = m_mapMethod == 0 ? strewn::KdeMethod::Adaptive : strewn::KdeMethod::Fixed;
	job.maps.cellM = number(m_values, "cell");
	if (m_mapMethod == 1)
		job.maps.fixedSigmaM = number(m_values, "sigma");
	job.maps.offline = m_offline;
	job.maps.basemap = strewn::Basemap::EsriImagery;
	if (m_environmentMode == 0) {
		strewn::AutomaticEnvironment environment;
		environment.offline = m_offline;
		job.environment = environment;
	} else if (m_environmentMode == 1)
		job.environment = strewn::EnvironmentFile{std::filesystem::path(m_environmentPath.toStdWString())};
	else
		job.environment = strewn::IdealizedEnvironment{number(m_values, "groundHeight")};
	strewn::validate(job.physics, job.integration);
	strewn::validate(mc);
	strewn::validate(job.maps);
	return job;
}

QString AtmosphereViewModel::sourceError() const
{
	if (m_sourceMode == 0) {
		if (m_loading)
			return tr("Reading and validating entry states...");
		if (!m_archiveError.isEmpty())
			return m_archiveError;
		if (m_archivePath.isEmpty())
			return tr("Select an Earth entry archive at 100 km.");
		if (!QFileInfo(m_archivePath).isFile())
			return tr("The selected archive is unavailable.");
		if (!m_archive)
			return tr("Select the entry archive again.");
		if (m_selectedEvent < 0)
			return tr("Select one entry event.");
		const auto& event = m_archive->events[m_selectedEvent];
		if (m_environmentMode == 0 && event.endMjdTdb - event.startMjdTdb > 1.0 / 24)
			return tr(
			    "Automatic weather supports entry events up to one hour. Use a saved environment or an explicitly "
			    "idealized "
			    "atmosphere for this event.");
		return {};
	}
	const auto text = m_values.value("epoch").toString();
	try {
		atmosphere::parseUtc(text);
		const auto latitude = number(m_values, "latitude"), longitude = number(m_values, "longitude");
		const auto speed = number(m_values, "speed"), heading = number(m_values, "heading"),
		           slope = number(m_values, "slope");
		if (std::abs(latitude) > 90 || std::abs(longitude) > 180)
			return tr("Latitude or longitude is outside its range.");
		if (speed <= 0 || heading < 0 || heading >= 360 || slope < -90 || slope >= 0)
			return tr("Require speed > 0, heading in [0, 360) and an inbound slope in [-90, 0) degrees.");
	} catch (const std::exception& e) {
		return QString::fromUtf8(e.what());
	}
	return {};
}

QString AtmosphereViewModel::formError() const
{
	try {
		parameterOptions();
	} catch (const std::exception& e) {
		return QString::fromUtf8(e.what());
	}
	if (m_environmentMode == 1 && !QFileInfo(m_environmentPath).isFile())
		return tr("Select a saved environment profile.");
	return {};
}
