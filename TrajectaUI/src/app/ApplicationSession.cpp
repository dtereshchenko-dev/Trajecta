#include "ApplicationSession.h"

#include <QDesktopServices>
#include <QFileInfo>
#include <QSettings>
#include <algorithm>
#include <cmath>
#include <optional>

namespace {
	const QStringList themeNames{"graphite", "observatory", "slate", "nocturne", "redlight", "dim", "light"};

	QString readableInputPath(const QUrl& url)
	{
		if (!url.isLocalFile() && !url.scheme().isEmpty())
			return {};
		const QFileInfo file(url.isLocalFile() ? url.toLocalFile() : url.toString());
		return file.isFile() && file.isReadable() ? file.absoluteFilePath() : QString();
	}
}  // namespace

ApplicationSession::ApplicationSession(FitSettingsAggregator& settings, QObject* parent)
    : QObject(parent), m_engine(settings), m_calculationSettings(settings), m_dataSources(&settings)
{
	m_dataSources.refresh();
	m_theme = QSettings().value("appearance/theme", "graphite").toString();
	if (!themeNames.contains(m_theme))
		m_theme = "graphite";
	connect(this, &ApplicationSession::changed, this, [this] {
		m_solution.observations()->setEditingEnabled(!busy());
		m_risk.setBlocked(m_busy || m_ephemerides.busy() || m_atmosphere.busy() || m_approaches.busy() || m_mov.busy());
		m_ephemerides.setBlocked(m_busy || m_risk.busy() || m_atmosphere.busy() || m_approaches.busy() || m_mov.busy());
		m_atmosphere.setBlocked(m_busy || m_risk.busy() || m_ephemerides.busy() || m_approaches.busy() || m_mov.busy());
		m_approaches.setBlocked(m_busy || m_risk.busy() || m_ephemerides.busy() || m_atmosphere.busy() || m_mov.busy());
		m_mov.setBlocked(m_busy || m_risk.busy() || m_ephemerides.busy() || m_atmosphere.busy() || m_approaches.busy());
	});
	connect(&m_risk, &RiskMonitorViewModel::changed, this, &ApplicationSession::changed);
	connect(&m_ephemerides, &EphemerisViewModel::changed, this, &ApplicationSession::changed);
	connect(&m_atmosphere, &AtmosphereViewModel::changed, this, &ApplicationSession::changed);
	connect(&m_approaches, &CloseApproachesViewModel::changed, this, &ApplicationSession::changed);
	connect(&m_mov, &MovViewModel::changed, this, &ApplicationSession::changed);
	connect(&m_engine, &OrbitEngine::progressStatusChanged, this, [this] {
		m_progress = m_cancelling ? tr("Cancelling...") : m_engine.progressStatus();
		emit changed();
	});
	connect(&m_engine, &OrbitEngine::errorMessageChanged, this, [this] {
		m_gaiaDownloading = false;
		m_error = m_cancelling ? QString() : m_engine.errorMessage();
		m_progress = m_cancelling ? tr("Cancelled") : tr("Calculation failed");
		m_busy = false;
		m_cancelling = false;
		m_preparedInput.reset();
		m_observationReview.clear();
		m_reviewingInput = false;
		if (!m_acceptedInputFile.isEmpty()) {
			m_inputFile = m_acceptedInputFile;
			m_loadRequested = false;
		}
		emit changed();
	});
	connect(&m_engine, &OrbitEngine::observationsReady, this,
	        [this](std::shared_ptr<orbital_core::ObservationSet> observations) {
		        m_preparedInput = std::move(observations);
		        m_observationReview.setObservations(*m_preparedInput);
		        m_busy = false;
		        m_cancelling = false;
		        if (m_observationReview.totalCount() > 0) {
			        m_reviewingInput = true;
			        m_progress = tr("Review observations before solving");
			        emit changed();
			        emit observationReviewReady();
		        } else {
			        run();
		        }
	        });
	connect(&m_engine, &OrbitEngine::solutionReady, this,
	        [this](const orbital_core::FitPipelineResult& result, const QString& directory) {
		        m_gaiaDownloading = false;
		        const auto moids = result.fit.orbit && !result.fit.cancelled
		                               ? orbital_core::hazard::computeMoidsToAllPlanets(*result.fit.orbit,
		                                                                                result.fit.orbit->getEpoch())
		                               : orbital_core::hazard::PlanetMoidTable{};
		        if (m_solution.publish(result, directory, moids)) {
			        const double magnitude = m_solution.absoluteMagnitude();
			        m_ephemerides.setSolution(
			            result, m_solution.objectTitle(),
			            std::isfinite(magnitude) ? std::optional<double>(magnitude) : std::nullopt,
			            m_solution.cometPhotometry() ? &m_solution.cometPhotometrySolution() : nullptr);
			        m_acceptedInputFile = m_inputFile;
			        m_loadRequested = false;
			        m_progress = tr("Result saved to working folder");
			        m_error.clear();
		        } else {
			        m_error = tr("The calculation returned no solution");
		        }
		        m_busy = false;
		        m_cancelling = false;
		        m_preparedInput.reset();
		        m_observationReview.clear();
		        emit changed();
	        });
}

QString ApplicationSession::inputName() const { return QFileInfo(m_inputFile).fileName(); }

bool ApplicationSession::downloadGaiaFpr()
{
	if (busy() || !m_solution.available())
		return false;
	const auto& snapshot = m_solution.snapshot();
	if (!snapshot || snapshot->observations.targetNaifId <= 2000000 || snapshot->observations.targetNaifId > 2999999) {
		m_error = tr("Gaia FPR download requires a numbered asteroid.");
	} else if (!m_engine.currentSettings().useWeights) {
		m_error = tr("Enable observation weighting before adding Gaia FPR.");
	} else {
		const double magnitude = m_solution.absoluteMagnitude();
		const double estimate =
		    std::isfinite(magnitude) ? 1329.0 * std::pow(10.0, -magnitude / 5.0) / (2.0 * std::sqrt(0.125)) : 0.0;
		const auto fallbackRadius =
		    std::isfinite(estimate) && estimate > 0.0 ? std::optional<double>(estimate) : std::nullopt;
		m_error.clear();
		m_busy = true;
		m_cancelling = false;
		m_gaiaDownloading = true;
		m_progress = tr("Downloading Gaia FPR observations...");
		m_engine.downloadGaiaFpr(snapshot->observations, m_solution.observations()->overrides(), fallbackRadius);
		emit changed();
		emit runStarted();
		return true;
	}
	emit changed();
	return false;
}

void ApplicationSession::setTheme(const QString& value)
{
	if (m_theme == value || !themeNames.contains(value))
		return;
	m_theme = value;
	QSettings().setValue("appearance/theme", value);
	emit themeChanged();
}

bool ApplicationSession::openInput(const QUrl& url)
{
	if (busy())
		return false;
	const auto path = readableInputPath(url);
	if (path.isEmpty()) {
		m_error = tr("Cannot read the observations file");
		emit changed();
		return false;
	}
	m_inputFile = path;
	m_loadRequested = true;
	m_preparedInput.reset();
	m_observationReview.clear();
	m_reviewingInput = false;
	m_recentFiles.addFile(QUrl::fromLocalFile(path));
	run();
	return true;
}

bool ApplicationSession::openRecent(int index)
{
	const auto path = m_recentFiles.filePath(index);
	if (path.isEmpty())
		return false;
	return openInput(QUrl::fromLocalFile(path));
}

bool ApplicationSession::canDropInput(const QList<QUrl>& urls) const
{
	return !busy() && urls.size() == 1 && urls.front().isLocalFile() && !readableInputPath(urls.front()).isEmpty();
}

void ApplicationSession::run()
{
	if (!canRun() || m_reviewingInput)
		return;
	m_busy = true;
	m_cancelling = false;
	m_error.clear();
	const bool reuseObservations = !m_loadRequested && m_solution.available();
	m_progress = reuseObservations ? tr("Recalculating solution...") : tr("Loading observations...");
	emit changed();
	emit runStarted();
	if (reuseObservations)
		m_engine.recalculate(m_solution.snapshot()->observations, m_solution.observations()->overrides());
	else if (m_preparedInput)
		m_engine.recalculate(*m_preparedInput, m_observationReview.overrides());
	else
		m_engine.prepareFile(m_inputFile);
}

void ApplicationSession::acceptObservationReview()
{
	if (!m_reviewingInput || !m_preparedInput || busy())
		return;
	m_reviewingInput = false;
	run();
}

void ApplicationSession::cancelObservationReview()
{
	if (!m_reviewingInput)
		return;
	m_reviewingInput = false;
	m_preparedInput.reset();
	m_observationReview.clear();
	m_inputFile = m_acceptedInputFile;
	m_loadRequested = m_acceptedInputFile.isEmpty();
	m_progress.clear();
	emit changed();
}

void ApplicationSession::cancel()
{
	if (m_gaiaDownloading) {
		m_cancelling = true;
		m_progress = tr("Cancelling...");
		m_engine.cancel();
		emit changed();
		return;
	}
	if (m_mov.busy()) {
		m_mov.cancel();
		return;
	}
	if (m_approaches.busy()) {
		m_approaches.cancel();
		return;
	}
	if (m_atmosphere.busy()) {
		m_atmosphere.cancel();
		return;
	}
	if (m_risk.busy()) {
		m_risk.cancel();
		return;
	}
	if (m_ephemerides.busy()) {
		m_ephemerides.cancel();
		return;
	}
	if (!m_busy || m_cancelling)
		return;
	m_cancelling = true;
	m_progress = tr("Cancelling...");
	m_engine.cancel();
	emit changed();
}

void ApplicationSession::openDirectory()
{
	if (m_solution.directory().isEmpty())
		return;
	if (!QDesktopServices::openUrl(QUrl::fromLocalFile(m_solution.directory()))) {
		m_error = tr("Cannot open the working folder");
		emit changed();
	}
}

void ApplicationSession::dismissError()
{
	m_error.clear();
	emit changed();
}

bool ApplicationSession::applySettings(bool recalculate)
{
	if (busy() || (recalculate && !canRun()) || !m_calculationSettings.apply())
		return false;
	if (recalculate)
		run();
	return true;
}
