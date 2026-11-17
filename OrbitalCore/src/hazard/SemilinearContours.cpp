#include "SemilinearContours.h"
#include "BoundaryScan.h"
#include "uncertainty/LinearizedCovariance.h"
#include "orbital_core/Constants.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include <Eigen/SVD>
#include <cmath>

namespace orbital_core::hazard {
	SemilinearContours computeSemilinearContours(const std::array<double, 6>& state, double epoch,
	                                             const NonGravModel& nonGrav, const FitParameterCovariance& covariance,
	                                             const std::optional<BoundaryCrossing>& nominal, Planet target,
	                                             double start, double end, double altitude, unsigned workers,
	                                             const NBodyPropagator& propagator, IEphemerisSource& ephemeris,
	                                             std::stop_token cancel)
	{
		SemilinearContours result;
		result.status = "unavailable";
		if (!nominal) {
			result.reason = "The nominal orbit does not intersect the selected boundary in the interval.";
			return result;
		}
		const int dimension = 6 + nonGrav.fitParameterCount();
		const auto factor = uncertainty::covarianceFactor(covariance, dimension);
		if (!factor) {
			result.reason = "A positive definite full fit-epoch covariance is required.";
			return result;
		}
		const Vec3D center = icrfToBodyFixedMjd(target, nominal->relativePositionKm, nominal->mjdTdb, &ephemeris);
		const Vec3D normal = center.normalized();
		if (!std::isfinite(center.dot(center)) || normal.isNull()) {
			result.reason = "The nominal body-fixed position is unavailable.";
			return result;
		}
		const Vec3D axis = std::abs(normal.z) < 0.9 ? Vec3D(0, 0, 1) : Vec3D(1, 0, 0);
		const Vec3D east = axis.cross(normal).normalized(), north = normal.cross(east);
		const auto trace = [&](const Eigen::VectorXd& offset) {
			auto perturbed = state;
			auto ng = nonGrav;
			for (int i = 0; i < 6; ++i)
				perturbed[i] += offset(i);
			for (int i = 0; i < ng.fitParameterCount(); ++i)
				ng.setFitParameterValue(i, ng.fitParameterValue(i) + offset(6 + i));
			auto p = propagator;
			return boundary::search(perturbed, epoch, ng, target, start, end, altitude, p, ephemeris, cancel);
		};
		try {
			// Differentiate the surface map in whitened fit parameters, then select
			// its two-dimensional covariance regression subspace with an SVD.
			const auto differentiate = [&](double step) {
				Eigen::MatrixXd jacobian(2, dimension);
				for (int i = 0; i < dimension; ++i) {
					const Eigen::VectorXd delta = step * factor->col(i);
					const auto plus = trace(delta), minus = trace(-delta);
					if (!plus.entry || !minus.entry)
						throw std::runtime_error(
						    "A perturbed orbit missed the boundary while differentiating the surface map");
					const Vec3D a =
					    icrfToBodyFixedMjd(target, plus.entry->relativePositionKm, plus.entry->mjdTdb, &ephemeris);
					const Vec3D b =
					    icrfToBodyFixedMjd(target, minus.entry->relativePositionKm, minus.entry->mjdTdb, &ephemeris);
					jacobian(0, i) = (a - b).dot(east) / (2 * step);
					jacobian(1, i) = (a - b).dot(north) / (2 * step);
				}
				return jacobian;
			};
			const Eigen::MatrixXd coarse = differentiate(0.05);
			const Eigen::MatrixXd derivative = differentiate(0.025);
			if (!derivative.allFinite())
				throw std::runtime_error("Non-finite surface-map derivative");
			if ((coarse - derivative).norm() > 0.1 * derivative.norm())
				throw std::runtime_error(
				    "Surface-map derivatives did not converge when the difference step was halved");
			const Eigen::JacobiSVD<Eigen::MatrixXd> svd(derivative, Eigen::ComputeThinV);
			if (svd.info() != Eigen::Success || svd.singularValues()(0) <= 0 ||
			    svd.singularValues()(1) <= svd.singularValues()(0) * 1e-10)
				throw std::runtime_error("The local surface map is rank deficient");
			const Eigen::MatrixXd basis = *factor * svd.matrixV().leftCols(2);
			constexpr int samples = 72;
			for (const int sigma : {1, 3, 5}) {
				SemilinearContour contour;
				contour.sigma = sigma;
				contour.points.resize(samples);
				std::vector<int> outcomes(samples, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
				for (int i = 0; i < samples; ++i) {
					if (cancel.stop_requested()) {
						outcomes[i] = 2;
						continue;
					}
					try {
						const double angle = Constants::TWOPI * i / samples;
						const Eigen::Vector2d direction(std::cos(angle), std::sin(angle));
						const auto hit = trace(sigma * basis * direction);
						contour.points[i] = hit.entry;
						if (!hit.entry)
							outcomes[i] = 1;
					} catch (...) {
						outcomes[i] = 2;
					}
				}
				for (const auto outcome : outcomes) {
					contour.missed += outcome == 1;
					contour.failed += outcome == 2;
				}
				result.contours.push_back(std::move(contour));
				if (cancel.stop_requested())
					throw std::runtime_error("Contour calculation cancelled");
			}
			result.status = "complete";
		} catch (const std::exception& e) {
			result.reason = e.what();
			result.status = cancel.stop_requested() ? "cancelled" : "unavailable";
			result.contours.clear();
		}
		return result;
	}
}  // namespace orbital_core::hazard
