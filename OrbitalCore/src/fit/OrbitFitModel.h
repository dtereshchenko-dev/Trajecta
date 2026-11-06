#pragma once

#include "fit/OrbitFitConfig.h"
#include "orbital_core/Orbit.h"
#include "NBodyPropagator.h"
#include "obs/PreparedObservation.h"
#include "obs/AstrometricResidual.h"
#include "orbital_core/nongrav/NonGravProfile.h"

namespace orbital_core::fit {

	// Owns fit-ready observations (see FitObservations.h) and the dynamical
	// model, not solver state.  It never prepares observations itself, so fits of
	// several subsets share one preparation.  Parameters use canonical velocities;
	// propagation uses AU/day and MJD TDB.
	// N-body residuals, variational equations and uncertainty transport share the
	// encounter controller (npoint=100); pure Kepler evaluation remains unchanged.
	class OrbitFitModel {
	   public:
		OrbitFitModel(Orbit_p initial, std::vector<obs::PreparedObservation> observations,
		              const OrbitOptimizerConfig& config);
		OrbitFitModel(const OrbitFitModel&) = delete;
		OrbitFitModel& operator=(const OrbitFitModel&) = delete;

		struct ResidualsAndStateJac {
			std::vector<double> residuals;
			Eigen::MatrixXd stateJ;
			std::vector<obs::AstrometricResidualTransform> transforms;
		};

		[[nodiscard]] const auto& observations() const noexcept { return observations_; }

		[[nodiscard]] int residualOffset(std::size_t i) const { return residualOffsets_.at(i); }

		[[nodiscard]] int numResiduals() const noexcept { return residualCount_; }

		[[nodiscard]] double opticalRms(const std::vector<double>& residuals) const;

		[[nodiscard]] bool usesWeights() const noexcept { return useWeights_; }

		[[nodiscard]] bool usesNBody() const noexcept { return useNBody_; }

		[[nodiscard]] const Orbit& initialOrbit() const noexcept { return *initialOrbit_; }

		[[nodiscard]] NBodyPropagator& propagator() const noexcept { return nbody_; }

		[[nodiscard]] const NonGravProfile& nonGravProfile() const noexcept { return nonGravProfile_; }

		[[nodiscard]] auto& differenceSteps() noexcept { return deltaVals_; }

		[[nodiscard]] const auto& differenceSteps() const noexcept { return deltaVals_; }

		[[nodiscard]] std::vector<double> computeResiduals(const std::vector<double>& params) const;
		[[nodiscard]] ResidualsAndStateJac computeResidualsAndStateJac(const std::vector<double>& params) const;
		[[nodiscard]] std::vector<Vec3D> computeLightTimePositions(const std::vector<double>& params,
		                                                           const NonGravModel& nonGrav) const;
		[[nodiscard]] Eigen::MatrixXd computeJacobian(const std::vector<double>& params) const;
		[[nodiscard]] double computeCost(const std::vector<double>& params) const;
		[[nodiscard]] double computeCostWeighted(const std::vector<double>& params) const;
		[[nodiscard]] Orbit_p orbitFromStateVector(const std::vector<double>& params) const;
		[[nodiscard]] static std::vector<double> stateVectorFromOrbit(const Orbit& orbit);

	   private:
		std::vector<NBodyPropagator::StateVel> receptionStates(const std::vector<double>& params,
		                                                       const NonGravModel& nonGrav) const;
		void configureDynamics(const OrbitOptimizerConfig& config);
		std::vector<double> computeResiduals(const std::vector<double>& params, const NonGravModel& nonGrav) const;
		Vec3D toTopocentric(const Vec3D& position, const obs::PreparedObservation& observation) const;
		static std::pair<double, double> toRADEC(const Vec3D& position) noexcept;

		Orbit_p initialOrbit_;
		std::vector<obs::PreparedObservation> observations_;
		std::vector<int> residualOffsets_;
		int residualCount_{};
		std::vector<double> deltaVals_{1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8};
		double cosObliquity_;
		double sinObliquity_;
		bool useWeights_;
		bool useNBody_;
		mutable NBodyPropagator nbody_;
		NonGravProfile nonGravProfile_{nbody_.nonGrav()};
	};

}  // namespace orbital_core::fit
