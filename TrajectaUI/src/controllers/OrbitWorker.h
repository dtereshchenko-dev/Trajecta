#pragma once

#include <QObject>
#include <QByteArray>
#include <cstdint>
#include <filesystem>
#include <variant>

#include "NBodyPropagator.h"
#include <stop_token>
#include "orbital_core/FitPipeline.h"

struct OptimizationSettings {
	orbital_core::NonGravModelKind nonGravModel{orbital_core::NonGravModelKind::MarsdenH2O};

	bool pertEnabled{false};
	bool useSb441{false};
	bool nonGravEnabled{false};
	bool debiasEnabled{false};
	bool useWeights{false};
	bool observerSigmas{false};
	bool timingEnabled{false};
	bool sigmaClipEnabled{false};
	orbital_core::PhotometryOptions photometry;
	orbital_core::ForceModelSettings forceModel;

	int nonGravSolveMask{0};

	orbital_core::AstrometricWeightProfile weightProfile{orbital_core::AstrometricWeightProfile::Vfcc17};

	uint16_t pertPlanetMask{static_cast<uint16_t>(PLANET_BODY_MASK | MOON_BODY_MASK)};

	double integratorInitialStepDays{10.0};
	double sigmaClipThreshold{3.0};
};

struct OrbitObservationInput {
	orbital_core::ObservationSet observations;
	std::vector<orbital_core::ObservationOverride> overrides;
};

using OrbitFitInput = std::variant<QString, OrbitObservationInput>;

struct OrbitFitJob {
	OptimizationSettings settings;
	OrbitFitInput input;
	QByteArray settingsSnapshot;
};

class OrbitWorker : public QObject {
	Q_OBJECT

   public:
	explicit OrbitWorker(QObject* parent = nullptr);
	~OrbitWorker() override = default;

	[[nodiscard]] const orbital_core::FitPipelineResult& result() const noexcept { return m_result; }

	[[nodiscard]] const std::filesystem::path& solutionDirectory() const noexcept { return m_solutionDirectory; }

	[[nodiscard]] static orbital_core::FitConfig makeFitConfig(const OptimizationSettings& settings);

   public slots:
	void process(const OrbitFitJob& job, std::stop_token cancel);

   signals:
	void optimizationFinished();
	void errorOccurred(const QString& message);
	void progressUpdate(const QString& status);

   private:
	orbital_core::FitPipelineResult m_result;
	std::filesystem::path m_solutionDirectory;
};
