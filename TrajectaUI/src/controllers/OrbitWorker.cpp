#include "controllers/OrbitWorker.h"

#include "orbital_core/ObservationSet.h"
#include "logger/Logger.h"
#include "orbital_core/Exporter.h"

#include <QCoreApplication>
#include <QFileInfo>
#include <QSaveFile>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {
	constexpr double kDefaultNBodyStepDays = 10.0;
	constexpr double kMinSaneNBodyStepDays = 1.0e-6;

	double saneNBodyStepDays(double stepDays)
	{
		if (!std::isfinite(stepDays) || stepDays < kMinSaneNBodyStepDays)
			return kDefaultNBodyStepDays;
		return stepDays;
	}

}  // namespace

OrbitWorker::OrbitWorker(QObject* parent) : QObject(parent) {}

orbital_core::FitConfig OrbitWorker::makeFitConfig(const OptimizationSettings& settings)
{
	namespace oc = orbital_core;
	oc::FitConfig cfg;

	cfg.useWeights = settings.useWeights;
	cfg.obsWeighting.useObserverSigmas = settings.observerSigmas;
	cfg.obsWeighting.timingSigmaDefaultSeconds = settings.timingEnabled ? oc::kDefaultTimingSigmaSeconds : -1.0;
	cfg.weightProfile = settings.weightProfile;

	cfg.sigmaClip.enabled = settings.sigmaClipEnabled;
	cfg.sigmaClip.thresholdNSigma = settings.sigmaClipThreshold;
	cfg.photometry = settings.photometry;

	cfg.nbody.enabled = settings.pertEnabled;
	cfg.nbody.stepSize = saneNBodyStepDays(settings.integratorInitialStepDays);
	cfg.nbody.planetMask = settings.pertPlanetMask;
	cfg.nbody.useSb441 = settings.useSb441;
	cfg.nbody.forceModel = settings.forceModel;
	if (cfg.nbody.stepSize != settings.integratorInitialStepDays) {
		LOG_INFO << "[NBody] invalid initial step " << settings.integratorInitialStepDays << "d; using "
		         << cfg.nbody.stepSize << "d";
	}

	cfg.nonGrav.enabled = settings.nonGravEnabled;
	cfg.nonGrav.model = settings.nonGravModel;
	cfg.nonGrav.solveMask = static_cast<std::uint8_t>(settings.nonGravSolveMask);

	return cfg;
}

void OrbitWorker::process(const OrbitFitJob& job, std::stop_token cancel)
{
	auto fail = [this](const QString& msg) { emit errorOccurred(msg); };

	try {
		m_solutionDirectory.clear();
		if (cancel.stop_requested())
			throw std::runtime_error("Cancelled");
		const auto cfg = makeFitConfig(job.settings);
		const orbital_core::ObservationPreparationOptions preparation{
		    .debias = job.settings.debiasEnabled,
		    .weightProfile = job.settings.weightProfile,
		    .runtimeDirectory = QCoreApplication::applicationDirPath().toStdString()};
		orbital_core::FitRunControl control{.cancel = cancel};
		control.progress = [this](const std::string& text) {
			const auto message = QString::fromStdString(text);
			emit progressUpdate(message);
		};
		orbital_core::ObservationSet observations;
		std::span<const orbital_core::ObservationOverride> overrides;
		if (const auto* file = std::get_if<QString>(&job.input)) {
			observations = orbital_core::loadObservations(std::filesystem::path(file->toStdWString()));
			orbital_core::prepareObservations(observations, preparation);
		} else {
			const auto& input = std::get<OrbitObservationInput>(job.input);
			observations = orbital_core::prepareObservationCopy(input.observations, preparation);
			overrides = input.overrides;
		}
		const auto& prepared = *observations.preparation;
		const auto message = prepared.debiasRequested
		                         ? QString("Debiased %1 / %2 observations (catalog bias)")
		                               .arg(prepared.correctedAll)
		                               .arg(observations.all.size())
		                         : QString("Debias off - using raw astrometry (%1 obs)").arg(observations.all.size());
		control.progress(message.toStdString());
		auto result = orbital_core::FitPipeline{cfg}.run({observations, std::nullopt, overrides}, control);
		if (result.fit.cancelled || cancel.stop_requested())
			throw std::runtime_error("Cancelled");
		if (!result.fit.orbit)
			throw std::runtime_error("The calculation returned no solution");
		if (result.fit.converged) {
			orbital_core::ExportContext exportContext;
			exportContext.objectName =
			    observations.objectName.empty() ? observations.source.stem().string() : observations.objectName;
			exportContext.fitResult = &result.fit;
			exportContext.nominalEarthEntry = result.nominalEarthEntry ? &*result.nominalEarthEntry : nullptr;
			exportContext.nominalEarthEntryError = result.nominalEarthEntryError;
			exportContext.observations = observations.all;
			exportContext.gaiaPhotocenterRadiusKm = observations.gaiaPhotocenterRadiusKm;
			exportContext.weightProfile = result.effectiveConfig.weightProfile;
			exportContext.useWeights = result.effectiveConfig.useWeights;
			exportContext.cometPhotometry = result.usesCometPhotometry();
			exportContext.photometry = result.effectiveConfig.photometry;

			const auto exportReport = orbital_core::Exporter::save(exportContext, {orbital_core::ExportProduct::All});
			QSaveFile settingsFile(QString::fromStdWString((exportReport.directory / "Trajecta.ini").wstring()));
			if (!settingsFile.open(QIODevice::WriteOnly) ||
			    settingsFile.write(job.settingsSnapshot) != job.settingsSnapshot.size() || !settingsFile.commit()) {
				throw std::runtime_error("Could not save the settings snapshot with the solution");
			}
			LOG_INFO << "Fit export written to " << exportReport.directory.string();
			m_solutionDirectory = exportReport.directory;
		}

		if (cancel.stop_requested())
			throw std::runtime_error("Cancelled");
		m_result = std::move(result);
		emit optimizationFinished();
	}

	catch (const std::exception& e) {
		fail(QString::fromStdString(e.what()));
	}

	catch (...) {
		fail(QStringLiteral("Unknown error"));
	}
}
