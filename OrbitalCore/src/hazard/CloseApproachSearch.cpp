#include "orbital_core/hazard/CloseApproachSearch.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/TimeScale.h"
#include "CloseApproachScan.h"
#include "OpikTargetPlane.h"
#include "TargetPlaneDistanceBounds.h"
#include "NBodyPropagator.h"
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard {
	namespace {
		constexpr double auKm = Constants::AU_IN_METERS / 1000.0;

		void checkCancelled(const FitRunControl& control)
		{
			if (control.cancel.stop_requested())
				throw std::runtime_error("Close-approach search cancelled");
		}

		std::optional<Eigen::MatrixXd> covarianceFactor(const FitResult& fit)
		{
			const int n = 6 + fit.finalNonGrav.fitParameterCount();
			Eigen::MatrixXd cov = fit.uncertainty.fitParameterCovariance;
			if (n == 6 && cov.size() == 0)
				cov = fit.uncertainty.cartesianStateCovariance;
			if (cov.rows() != n || cov.cols() != n || !cov.allFinite() || (cov.diagonal().array() <= 0.0).any())
				return std::nullopt;
			const Eigen::VectorXd scale = cov.diagonal().cwiseSqrt();
			const Eigen::MatrixXd correlation =
			    scale.cwiseInverse().asDiagonal() * cov * scale.cwiseInverse().asDiagonal();
			if (!correlation.isApprox(correlation.transpose(), 1e-10))
				return std::nullopt;
			const Eigen::LLT<Eigen::MatrixXd> llt(correlation);
			if (llt.info() != Eigen::Success)
				return std::nullopt;
			return scale.asDiagonal() * Eigen::MatrixXd(llt.matrixL());
		}

		void addUncertainty(NominalCloseApproach& row, const NBodyPropagator::PropagationResultFull& state,
		                    const Eigen::MatrixXd& factor, const NBodyPropagator& prop, const NonGravModel& ng,
		                    IEphemerisSource& eph)
		{
			const auto& event = row.encounter;
			const double mjd = event.jd - 2400000.5;
			const auto planet = eph.queryMjd(toEphemerisPlanet(event.planet), mjd);
			const Vec3D r = state.R - planet.R, v = state.V - planet.V;
			const double distance = r.magnitude();
			const double nominalDifference = (state.R - event.rTarget).magnitude();
			if (nominalDifference > std::max(1e-8, event.distAu * 1e-5)) {
				row.uncertaintyNote = "Variational and nominal trajectories disagree";
				return;
			}
			constexpr double dt = 0.001;
			const Vec3D planetAcceleration = (eph.queryMjd(toEphemerisPlanet(event.planet), mjd + dt).V -
			                                  eph.queryMjd(toEphemerisPlanet(event.planet), mjd - dt).V) /
			                                 (2.0 * dt);
			const Vec3D acceleration = prop.accelerationAtMjd(state.R, state.V, mjd, ng) - planetAcceleration;
			const double denominator = v.dot(v) + r.dot(acceleration);
			if (!(distance > 0.0) || !std::isfinite(denominator) ||
			    denominator <= 1e-10 * (v.dot(v) + std::abs(r.dot(acceleration)))) {
				row.uncertaintyNote = "Degenerate encounter minimum";
				return;
			}
			Eigen::MatrixXd jacobian(6, factor.rows());
			jacobian.leftCols<6>() = state.Phi;
			for (int i = 0; i < ng.fitParameterCount(); ++i)
				jacobian.col(6 + i) = state.SensA.col(ng.sensitivityColumnForFitParameter(i));
			const Eigen::MatrixXd weighted = jacobian * factor;
			const Eigen::RowVector3d radial(r.x, r.y, r.z), velocity(v.x, v.y, v.z);
			const Eigen::RowVectorXd timeGradient =
			    -(velocity * weighted.topRows<3>() + radial * weighted.bottomRows<3>()) / denominator;
			const double threeSigmaTime = 3.0 * timeGradient.norm();
			if (!std::isfinite(threeSigmaTime)) {
				row.uncertaintyNote = "Non-finite propagated covariance";
				return;
			}
			row.timeUncertaintyDays = threeSigmaTime;

			const Eigen::Matrix3d toEquatorial =
			    Eigen::AngleAxisd(Constants::EARTH_INCL, Eigen::Vector3d::UnitX()).toRotationMatrix();
			const Eigen::Vector3d rEq = toEquatorial * radial.transpose() * auKm;
			const Eigen::Vector3d vEq = toEquatorial * velocity.transpose() * (auKm / 86400.0);
			const auto plane = computeOpikTargetPlaneDifferential({rEq.x(), rEq.y(), rEq.z()},
			                                                      {vEq.x(), vEq.y(), vEq.z()}, event.planet);
			if (!plane.valid) {
				row.uncertaintyNote = "B-plane distance bounds unavailable; local linear time uncertainty only";
				return;
			}
			Eigen::MatrixXd atMinimum(6, weighted.cols());
			atMinimum.topRows<3>() =
			    toEquatorial * (weighted.topRows<3>() + velocity.transpose() * timeGradient) * auKm;
			atMinimum.bottomRows<3>() =
			    toEquatorial *
			    (weighted.bottomRows<3>() +
			     Eigen::Vector3d(acceleration.x, acceleration.y, acceleration.z) * timeGradient) *
			    (auKm / 86400.0);
			const Eigen::MatrixXd projected = plane.dCoordsDGeoState * atMinimum;
			if (!projected.allFinite()) {
				row.uncertaintyNote = "Non-finite B-plane covariance; local linear time uncertainty only";
				return;
			}
			const Eigen::JacobiSVD<Eigen::MatrixXd> svd(projected, Eigen::ComputeThinU);
			const Eigen::Vector2d semiAxes = 3.0 * svd.singularValues();
			const auto bounds =
			    detail::ellipseDistanceBounds(svd.matrixU().transpose() * plane.point.coordsKm, semiAxes);
			if (!bounds) {
				row.uncertaintyNote = "Invalid B-plane ellipse; local linear time uncertainty only";
				return;
			}
			const double alpha = planetMuKm3PerS2(event.planet) / (plane.point.vInfKmS * plane.point.vInfKmS);
			row.minimumDistanceAu = detail::pericentreFromImpactParameter(bounds->minimum, alpha) / auKm;
			row.maximumDistanceAu = detail::pericentreFromImpactParameter(bounds->maximum, alpha) / auKm;
			row.uncertaintyNote = "Linearized 3-sigma B-plane ellipse; nominal-energy focusing";
			if (semiAxes.x() >= plane.point.coordsKm.norm() || threeSigmaTime > 10.0)
				row.uncertaintyNote += "; large uncertainty, nonlinear verification required";
		}
	}  // namespace

	double defaultCloseApproachLimitAu(Planet planet) noexcept
	{
		switch (planet) {
			case Planet::Venus:
				return 0.15;
			case Planet::Earth:
				return 0.5;
			case Planet::Jupiter:
			case Planet::Saturn:
			case Planet::Uranus:
			case Planet::Neptune:
				return 2.0;
			case Planet::Moon:
				return 0.01;
			default:
				return 0.1;
		}
	}

	CloseApproachSearchResult searchCloseApproaches(const FitPipelineResult& solution,
	                                                const CloseApproachSearchRequest& request,
	                                                const FitRunControl& control)
	{
		CloseApproachSearchResult result;
		try {
			checkCancelled(control);
			const auto& fit = solution.fit;
			const auto& forces = solution.effectiveConfig.nbody;
			if (!fit.fitEpochOrbit || fit.cancelled || !fit.converged || !forces.enabled)
				throw std::invalid_argument("Close approaches require a converged N-body solution");
			const double epoch = fit.fitEpochOrbit->epochMjd();
			if (!std::isfinite(epoch) || !std::isfinite(request.startMjdTdb) || !std::isfinite(request.endMjdTdb) ||
			    request.endMjdTdb <= request.startMjdTdb ||
			    std::max(std::abs(request.startMjdTdb - epoch), std::abs(request.endMjdTdb - epoch)) > 365250.0 ||
			    (request.maximumDistanceAu &&
			     (!std::isfinite(*request.maximumDistanceAu) || *request.maximumDistanceAu <= 0.0)) ||
			    !std::isfinite(forces.stepSize) || forces.stepSize <= 0.0)
				throw std::invalid_argument(
				    "Invalid close-approach interval or distance (maximum 1000 years from fit epoch)");
			auto* eph = getDefaultEphemerisSource();
			if (!eph || !eph->isAvailable())
				throw std::runtime_error("Planetary ephemerides are unavailable");
			if (forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
				throw std::runtime_error("SB441 is unavailable");
			NBodyPropagator prop;
			prop.setStepSize(forces.stepSize);
			prop.setPlanetMask(forces.planetMask);
			prop.setUseSb441(forces.useSb441);
			prop.setForceModel(forces.forceModel);
			prop.setTargetNaifId(forces.targetNaifId);
			prop.setInactiveInteriorMassPolicy(forces.inactiveInteriorMassPolicy);
			prop = prop.forEncounterMonitoring(*eph);
			prop.nonGrav() = fit.finalNonGrav;
			const Vec3D r = fit.fitEpochOrbit->getR();
			const Vec3D v = fit.fitEpochOrbit->getV() * Constants::GAUSS_K;
			double firstFutureImpact = std::numeric_limits<double>::infinity();
			double firstPastImpact = -std::numeric_limits<double>::infinity();
			for (Planet body : kAllPlanets) {
				checkCancelled(control);
				if (control.progress)
					control.progress("Searching " + std::string(longLabel(body)) + "...");
				CloseApproachConfig config;
				config.planet = body;
				config.horizonYears = std::max(1.0, request.endMjdTdb + 1.0 - epoch) / 365.25;
				config.lookbackYears = std::max(0.0, epoch - request.startMjdTdb + 1.0) / 365.25;
				config.thresholdAu = request.maximumDistanceAu.value_or(defaultCloseApproachLimitAu(body));
				config.encounterRadiusAu = config.thresholdAu;
				config.sampleStepDay = 0.25;
				const auto events =
				    findCloseApproachesChecked(r, v, epoch + 2400000.5, prop, fit.finalNonGrav, config, eph, control);
				for (const auto& event : events) {
					const double mjd = event.jd - 2400000.5;
					if (event.distAu * auKm < planetRadiusKm(body)) {
						if (mjd >= epoch)
							firstFutureImpact = std::min(firstFutureImpact, mjd);
						else
							firstPastImpact = std::max(firstPastImpact, mjd);
						result.stoppedAtSurface = true;
					}
					if (mjd < request.startMjdTdb || mjd > request.endMjdTdb)
						continue;
					NominalCloseApproach row;
					row.encounter = event;
					row.mjdTdb = mjd;
					row.relativeSpeedKmS = (event.vTarget - event.vPlanet).magnitude() * auKm / 86400.0;
					row.uncertaintyNote = "Fit-epoch covariance unavailable or invalid";
					result.events.push_back(std::move(row));
				}
			}
			std::erase_if(result.events, [&](const auto& row) {
				const double mjd = row.encounter.jd - 2400000.5;
				return mjd > firstFutureImpact || mjd < firstPastImpact;
			});
			std::sort(result.events.begin(), result.events.end(),
			          [](const auto& a, const auto& b) { return a.encounter.jd < b.encounter.jd; });
			const auto factor = covarianceFactor(fit);
			if (factor) {
				for (int direction : {-1, 1}) {
					checkCancelled(control);
					std::vector<std::size_t> indices;
					for (std::size_t i = 0; i < result.events.size(); ++i) {
						auto& row = result.events[i];
						if (row.encounter.distAu * auKm <= planetRadiusKm(row.encounter.planet)) {
							row.uncertaintyNote = "Interior osculating pericentre; uncertainty unavailable";
							continue;
						}
						if ((row.encounter.jd - 2400000.5 < epoch) == (direction < 0))
							indices.push_back(i);
					}
					if (direction < 0)
						std::reverse(indices.begin(), indices.end());
					if (indices.empty())
						continue;
					if (control.progress)
						control.progress("Propagating encounter covariance...");
					std::vector<double> times;
					for (auto i : indices)
						times.push_back(result.events[i].encounter.jd - 2400000.5);
					const auto states = prop.propagateAllWithStmAndSensMjd(r, v, epoch, times, fit.finalNonGrav,
					                                                       nullptr, control.cancel);
					checkCancelled(control);
					if (states.size() != indices.size())
						throw std::runtime_error("Incomplete covariance propagation");
					for (std::size_t i = 0; i < indices.size(); ++i)
						addUncertainty(result.events[indices[i]], states[i], *factor, prop, fit.finalNonGrav, *eph);
				}
			}
			checkCancelled(control);
			return result;
		} catch (...) {
			if (!control.cancel.stop_requested())
				throw;
			return {{}, true, false};
		}
	}
}  // namespace orbital_core::hazard
