#pragma once

#include "orbital_core/FitResult.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/BPlane.h"
#include "orbital_core/hazard/CloseApproach.h"

#include <Eigen/Core>

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

namespace orbital_core::hazard {

	struct ORBITAL_CORE_API ImpactRiskRmsScaling {
		using Objective = std::function<double(const Eigen::VectorXd&)>;

		struct ObservationLinearization {
			std::vector<double> residuals;
			Eigen::MatrixXd jacobian;
			double sumSquares{0.0};
			bool valid{false};
		};

		using Linearization = std::function<ObservationLinearization(const Eigen::VectorXd&)>;

		struct GammaNormalization {
			std::once_flag once;
			double gamma{1.0};
			double qNom{std::numeric_limits<double>::quiet_NaN()};
			double sampleMeanDelta{0.0};
			double sampleMeanS2{0.0};
			double minDeltaQ{0.0};
			double maxDeltaQ{0.0};
			int samples{0};
			bool valid{false};
		};

		Objective objective;
		Linearization linearization;
		double nominalChiSquare{std::numeric_limits<double>::quiet_NaN()};
		double objectiveScale{1.0};
		int numSamples{1000};
		std::shared_ptr<GammaNormalization> gammaCache;

		[[nodiscard]] bool enabled() const noexcept { return static_cast<bool>(objective) && numSamples > 0; }
	};

	struct ORBITAL_CORE_API ImpactRiskConfig {
		int maxIterations{30};
		double deltaXTolerance{1.0e-12};
		double alphaInitial{10.0};
		ImpactRiskRmsScaling rmsScaling{};
	};

	using ImpactRiskStateVec = Eigen::Matrix<double, 6, 1, Eigen::DontAlign>;
	using ImpactRiskStateMatrix = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;
	using ImpactRiskBPlaneVec = Eigen::Matrix<double, 2, 1, Eigen::DontAlign>;

	enum class ImpactRiskStatus : int {
		NotComputed = 0,
		CovarianceInvalid,
		EphemerisUnavailable,
		ReferencePropagationFailed,
		EncounterPlaneInvalid,
		NominalCovarianceInvalid,
		IterationPropagationFailed,
		EncounterMapInvalid,
		LinearSolveFailed,
		ViCovarianceInvalid,
		BeyondLinearRegime,
		Completed,
	};

	enum class ImpactProbabilitySource : int {
		None = 0,
		LineOfVariation,
	};

	struct ORBITAL_CORE_API ImpactRiskResult {
		ImpactRiskStatus status{ImpactRiskStatus::NotComputed};
		ImpactProbabilitySource probabilitySource{ImpactProbabilitySource::None};
		double probabilityImpact{0.0};
		double sigmaProbability{0.0};
		double sigmaLov{0.0};
		int numIterations{0};
		bool converged{false};
		double gamma{1.0};
		double vInfKmS{0.0};
		ImpactRiskStateVec xVI{ImpactRiskStateVec::Zero()};
		ImpactRiskStateMatrix sigmaVI{ImpactRiskStateMatrix::Zero()};
		ImpactRiskBPlaneVec bPlaneVIKm{ImpactRiskBPlaneVec::Zero()};
		BPlaneBasis basis{};
		std::vector<BPlanePoint> scatterKm{};
		BPlaneCovariance scatterCov{};
		double sigmaAlongMissKm{};
		double sigmaTransverseKm{};
		double sigmaTcaDays{std::numeric_limits<double>::quiet_NaN()};
		double minMissKm{};
		double medianMissKm{};
		int numImpactHits{};
		double bestMatchTcaJd{};
		double mahalanobisToNominal{0.0};
		bool beyondLinearRegime{false};
		double probabilityBeforeRmsScaling{0.0};
		double sigmaBeforeRmsScaling{0.0};
		double rmsScalingGamma{1.0};
		double rmsScalingCorrection{1.0};
		double rmsScalingDeltaChiSquare{0.0};
		double rmsScalingSigmaSquare{0.0};
		int rmsScalingSamples{0};
		bool rmsScalingApplied{false};
		bool computed{false};
		Eigen::VectorXd parameterVI;
		Eigen::VectorXd parameterStepPerSigmaVI;
		Eigen::MatrixXd parameterCovarianceVI;
	};

	struct ORBITAL_CORE_API ImpactRiskCluster {
		double startJd{};
		double endJd{};
		double meanJd{};
		double minDistAu{};
		double minDistLd{};
		int numEncounters{};
		int numSamples{};
		int firstEncounterIndex{-1};
		int lastEncounterIndex{-1};
		int representativeEncounterIndex{-1};
	};

	struct ORBITAL_CORE_API ImpactRiskRow {
		std::size_t clusterIndex{};
		std::size_t pathIndex{0};
		std::size_t pathCount{1};
		int representativeSampleIndex{-1};
		ImpactRiskCluster cluster{};
		CloseApproachEvent representativeEvent{};
		Eigen::VectorXd representativeParameters{};
		ImpactRiskResult result{};
	};

}  // namespace orbital_core::hazard
