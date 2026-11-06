#pragma once

#include "fit/OrbitFitProblem.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/FitResult.h"
#include <stop_token>

namespace orbital_core::fit {

	struct FitCovariance {
		Eigen::MatrixXd fullCovariance;
		Eigen::MatrixXd formalFullCovariance;
		double varianceFactor{0.0};
	};

	// Normal matrix at the fit epoch and derived elements at the output epoch.
	// The physical state covariance is never silently moved to the output epoch.
	class OrbitFitUncertainty {
	   public:
		OrbitFitUncertainty(const OrbitFitModel& model, const OrbitFitProblem& problem, double outputEpochTdb)
		    : model_(model), problem_(problem), outputEpochTdb_(outputEpochTdb)
		{
		}

		FitCovariance computeFitCovariance(const std::vector<double>& params, const std::vector<double>& rawResiduals,
		                                   bool log = true) const;

		struct ElementDerivatives {
			Eigen::Matrix<double, 6, Eigen::Dynamic> elemJ;
			Eigen::Matrix<double, 9, Eigen::Dynamic> qualityElemJ;
		};

		// Local element conversion is differentiated only in the six Cartesian
		// coordinates at the output epoch, then chained with STM/NG sensitivities.
		// No finite-difference parameter trial is propagated.
		ElementDerivatives computeElementDerivatives(const std::vector<double>& params, const Orbit& orbit,
		                                             std::stop_token cancel) const;

		// COT is differentiated analytically at the fit epoch, with canonical
		// velocity scaling and zero direct non-grav columns. Display transport
		// retains every fitted column. Quality rows are (q, 1/a, Tp), chained
		// through COT; no perturbed a or Q is evaluated. Failed transport or an
		// unavailable perihelion branch invalidates the block, never its epoch.
		struct CometaryDerivatives {
			Eigen::Matrix<double, 6, Eigen::Dynamic> cotJ;
			Eigen::Matrix<double, 3, Eigen::Dynamic> regularQualityJ;
			double epochJdTdb{std::numeric_limits<double>::quiet_NaN()};
			bool valid{false};
		};

		CometaryDerivatives computeCometaryDerivatives(const std::vector<double>& params, const Orbit& fitOrbit,
		                                               const Orbit& displayOrbit, std::stop_token cancel) const;
		// d(q,1/a,Tp)/d(COT). All entries are analytic except d(Tp)/de, which uses
		// a central difference with a representable, geometry-scaled step. Tp's
		// offset from the display epoch avoids JD cancellation. The nearest nominal
		// perihelion selects one local branch across aphelion; e is never clamped.
		using RegularQualityJacobian = Eigen::Matrix<double, 3, 6>;
		static std::optional<RegularQualityJacobian> regularQualityJacobianCot(const cometary::Vector6& nominalCot);
		// Projects full fitted covariance without modifying its Cartesian blocks.
		// Conventional covariance is derived analytically only where a is finite.
		static void populateCometaryUncertainty(Uncertainty& uncertainty, const CometaryDerivatives& derivatives,
		                                        const FitCovariance& covariance, const Orbit& displayOrbit);

		static uncertainty::DerivedSigmas computeDerivedSigmasFromQualityCovariance(
		    const Orbit& orbit, const uncertainty::ElementCov6x6& covariance, const Eigen::MatrixXd& qualityCovariance);
		static void populateOrbitElementSigmas(Orbit& orbit, const uncertainty::DerivedSigmas& sigmas) noexcept;
		static double orbitQualityU(const Orbit& orbit, const uncertainty::DerivedSigmas& sigmas) noexcept;
		static void logOrbitQualityU(const Orbit& orbit, const uncertainty::DerivedSigmas& sigmas);

	   private:
		// State and derivative at the output epoch in canonical velocity units.
		// Columns retain the fitted parameter order, including sparse A/DT masks.
		// Propagation uses variational equations, never perturbed NG trajectories.
		struct OutputTransport {
			std::vector<double> state;
			Eigen::Matrix<double, 6, Eigen::Dynamic> jacobian;
		};

		OutputTransport outputTransport(const std::vector<double>& params, double fitEpochTdb, double outputEpochTdb,
		                                std::stop_token cancel) const;
		// Unwrap elliptic samples to the nominal anomaly branch using each
		// sample's own period. An unbound sample cannot acquire an extra orbit.
		static double timeSincePerihelionOnNominalBranch(const cometary::Vector6& sampleCot,
		                                                 const cometary::Vector6& nominalCot);

		struct OrbitQualityMetrics {
			double u{std::numeric_limits<double>::quiet_NaN()};
			double runoff{std::numeric_limits<double>::quiet_NaN()};
			double termTp{0.0};
			double termInvA{0.0};
			bool valid{false};
			bool elliptical{false};
		};

		static OrbitQualityMetrics computeOrbitQualityMetrics(const Orbit& orbit,
		                                                      const uncertainty::DerivedSigmas& sigmas) noexcept;
		static std::vector<double> orbitalElements(const Orbit& orbit);
		static std::vector<double> qualityCovarianceElements(const Orbit& orbit);

		const OrbitFitModel& model_;
		const OrbitFitProblem& problem_;
		double outputEpochTdb_;
	};

}  // namespace orbital_core::fit
