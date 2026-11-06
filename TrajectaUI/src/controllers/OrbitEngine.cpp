#include "OrbitEngine.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/obs/GaiaFpr.h"
#include "orbital_core/obs/SbdbPhysical.h"

#include <QUrl>
#include <QStandardPaths>
#include <cstdint>
#include <filesystem>

OrbitEngine::OrbitEngine(FitSettingsAggregator& settings, QObject* parent)
    : QObject(parent), m_worker(new OrbitWorker), m_settings(settings)
{
	connect(&m_settings, &FitSettingsAggregator::anyChanged, this, [this] { m_settings.save(); });
	m_worker->moveToThread(&m_workerThread);
	connect(&m_workerThread, &QThread::finished, m_worker, &QObject::deleteLater);
	connect(m_worker, &OrbitWorker::optimizationFinished, this, [this] {
		m_fitRunning = false;
		emit solutionReady(m_worker->result(), QString::fromStdWString(m_worker->solutionDirectory().wstring()));
	});
	connect(m_worker, &OrbitWorker::errorOccurred, this, [this](const QString& message) {
		m_fitRunning = false;
		m_errorMessage = message;
		emit errorMessageChanged();
	});
	connect(m_worker, &OrbitWorker::progressUpdate, this, [this](const QString& message) {
		m_progressStatus = message;
		emit progressStatusChanged();
	});
	m_workerThread.start();
}

OrbitEngine::~OrbitEngine()
{
	disconnect(m_worker, nullptr, this, nullptr);
	m_fitCancel.request_stop();
	m_workerThread.quit();
	m_workerThread.wait();
}

void OrbitEngine::loadFile(const QString& path)
{
	if (m_fitRunning)
		return;
	submit(path.startsWith("file:///") ? QUrl(path).toLocalFile() : path);
}

void OrbitEngine::prepareFile(const QString& path)
{
	if (m_fitRunning)
		return;
	m_fitRunning = true;
	m_fitCancel = std::stop_source{};
	const auto cancel = m_fitCancel.get_token();
	QMetaObject::invokeMethod(
	    m_worker,
	    [this, path, cancel] {
		    try {
			    auto observations = std::make_shared<orbital_core::ObservationSet>(
			        orbital_core::loadObservations(std::filesystem::path(path.toStdWString())));
			    if (cancel.stop_requested())
				    throw std::runtime_error("Cancelled");
			    QMetaObject::invokeMethod(
			        this,
			        [this, observations = std::move(observations), cancel] {
				        m_fitRunning = false;
				        if (cancel.stop_requested()) {
					        m_errorMessage = QStringLiteral("Cancelled");
					        emit errorMessageChanged();
				        } else {
					        emit observationsReady(observations);
				        }
			        },
			        Qt::QueuedConnection);
		    } catch (const std::exception& error) {
			    const auto message = QString::fromUtf8(error.what());
			    QMetaObject::invokeMethod(
			        this,
			        [this, message] {
				        m_fitRunning = false;
				        m_errorMessage = message;
				        emit errorMessageChanged();
			        },
			        Qt::QueuedConnection);
		    }
	    },
	    Qt::QueuedConnection);
}

void OrbitEngine::recalculate(const orbital_core::ObservationSet& observations,
                              std::vector<orbital_core::ObservationOverride> overrides)
{
	if (m_fitRunning)
		return;
	submit(OrbitObservationInput{observations, std::move(overrides)});
}

void OrbitEngine::downloadGaiaFpr(orbital_core::ObservationSet observations,
                                  std::vector<orbital_core::ObservationOverride> overrides,
                                  std::optional<double> photocenterRadiusKm)
{
	if (m_fitRunning)
		return;
	m_fitRunning = true;
	m_fitCancel = std::stop_source{};
	const auto cancel = m_fitCancel.get_token();
	const auto cacheLocation = QStandardPaths::writableLocation(QStandardPaths::CacheLocation);
#ifdef _WIN32
	const std::filesystem::path cacheDirectory = cacheLocation.isEmpty()
	                                                 ? std::filesystem::path{}
	                                                 : std::filesystem::path(cacheLocation.toStdWString()) / "gaia";
#else
	const std::filesystem::path cacheDirectory =
	    cacheLocation.isEmpty() ? std::filesystem::path{} : std::filesystem::path(cacheLocation.toStdString()) / "gaia";
#endif
	QMetaObject::invokeMethod(
	    m_worker,
	    [this, observations = std::move(observations), overrides = std::move(overrides), photocenterRadiusKm,
	     cacheDirectory, cancel]() mutable {
		    try {
			    auto radius = photocenterRadiusKm;
			    auto radiusSource =
			        radius ? orbital_core::GaiaRadiusSource::FitEstimate : orbital_core::GaiaRadiusSource::None;
			    try {
				    if (const auto sbdbRadius =
				            orbital_core::obs::downloadSbdbRadiusKm(observations.targetNaifId - 2000000, cancel)) {
					    radius = sbdbRadius;
					    radiusSource = orbital_core::GaiaRadiusSource::Sbdb;
				    }
			    } catch (const std::exception&) {
				    if (cancel.stop_requested())
					    throw;
			    }
			    const auto gaia =
			        orbital_core::obs::downloadGaiaFpr(observations.targetNaifId - 2000000, cacheDirectory, cancel);
			    if (cancel.stop_requested())
				    throw std::runtime_error("Cancelled");
			    auto merged = orbital_core::withGaiaFprObservations(observations, gaia.observations, radius);
			    merged.gaiaRadiusSource = radiusSource;
			    const bool cached = gaia.source == orbital_core::obs::GaiaFprSource::Cache;
			    if (cached)
				    merged.format += " (cached)";
			    else if (!gaia.notice.empty())
				    merged.format += " (cache not saved)";
			    QMetaObject::invokeMethod(
			        this,
			        [this, merged = std::move(merged), overrides = std::move(overrides), cancel, cached]() mutable {
				        m_fitRunning = false;
				        if (cancel.stop_requested()) {
					        m_errorMessage = QStringLiteral("Cancelled");
					        emit errorMessageChanged();
					        return;
				        }
				        m_progressStatus =
				            cached ? tr("Solving with cached Gaia FPR...") : tr("Solving with Gaia FPR...");
				        emit progressStatusChanged();
				        submit(OrbitObservationInput{std::move(merged), std::move(overrides)});
			        },
			        Qt::QueuedConnection);
		    } catch (const std::exception& error) {
			    const auto message = QString::fromUtf8(error.what());
			    QMetaObject::invokeMethod(
			        this,
			        [this, message] {
				        m_fitRunning = false;
				        m_errorMessage = message;
				        emit errorMessageChanged();
			        },
			        Qt::QueuedConnection);
		    }
	    },
	    Qt::QueuedConnection);
}

void OrbitEngine::submit(OrbitFitInput input)
{
	m_fitRunning = true;
	m_fitCancel = std::stop_source{};
	OrbitFitJob job{currentSettings(), std::move(input)};
	QMetaObject::invokeMethod(
	    m_worker,
	    [worker = m_worker, job = std::move(job), cancel = m_fitCancel.get_token()] { worker->process(job, cancel); },
	    Qt::QueuedConnection);
}

OptimizationSettings OrbitEngine::currentSettings() const
{
	const auto* pert = m_settings.perturbations();
	const auto* ng = m_settings.nonGrav();
	const auto* w = m_settings.weighting();
	const auto* sc = m_settings.sigmaClip();
	const auto* photometry = m_settings.photometry();

	OptimizationSettings s;

	s.pertEnabled = pert->enabled() && orbital_core::planetaryEphemerisAvailable();
	s.integratorInitialStepDays = pert->integratorInitialStepDays();
	s.pertPlanetMask = static_cast<uint16_t>(pert->planetMask());
	s.useSb441 = pert->useSb441() && orbital_core::asteroidEphemerisAvailable();

	s.nonGravEnabled = ng->enabled();
	s.nonGravModel = static_cast<orbital_core::NonGravModelKind>(ng->model());
	s.nonGravSolveMask = ng->solveMask();

	s.useWeights = w->useWeights();
	s.observerSigmas = w->observerSigmas();
	s.timingEnabled = w->timingEnabled();
	s.debiasEnabled = w->debiasEnabled();
	const int weightProfile = w->weightProfile();
	s.weightProfile = static_cast<orbital_core::AstrometricWeightProfile>(
	    (weightProfile >= 0 && weightProfile <= 2) ? weightProfile : 1);

	s.sigmaClipEnabled = sc->enabled();
	s.sigmaClipThreshold = sc->threshold();
	s.photometry = photometry->options();

	return s;
}

void OrbitEngine::cancel() { m_fitCancel.request_stop(); }
