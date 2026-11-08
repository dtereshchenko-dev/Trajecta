#include "orbital_core/Ephemeris.h"

#include "EphemeridesContainer.h"
#include "TwoBodyPropagation.h"
#include "NBodyPropagator.h"
#include "obs/LightDeflection.h"
#include "obs/LightTime.h"
#include "obs/LightTimeSeries.h"
#include "obs/ObservatoryCatalog.h"
#include "obs/ObserverGeometry.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core {
	namespace {

		struct ObserverSite {
			bool ground{false};
			obs::ObserverData data;
		};

		bool finite(const Vec3D& v) noexcept { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

		ObserverSite resolveObserver(const EphemerisRequest& request)
		{
			ObserverSite site;
			if (request.center != EphemerisCenter::Station)
				return site;
			const auto& catalog = obs::observatoryCatalog();
			const auto station = catalog.find(request.station);
			if (station == catalog.end())
				throw std::invalid_argument("Unknown observatory code: " + request.station);
			const auto& record = station->second;
			if (record.code != "500" && record.cosPhi == 0.0 && record.sinPhi == 0.0)
				throw std::invalid_argument("Observatory " + request.station +
				                            " has no fixed geodetic position (spacecraft)");
			site.ground = true;
			site.data = obs::ObserverData{record.longitude, record.cosPhi, record.sinPhi, record.code};
			return site;
		}

		NBodyPropagator ephemerisPropagator(const EphemerisSolution& solution)
		{
			NBodyPropagator propagator;
			propagator.setStepSize(solution.forces.stepSize);
			propagator.setPlanetMask(solution.forces.enabled ? solution.forces.planetMask : 0);
			propagator.setUseSb441(solution.forces.useSb441);
			propagator.setTargetNaifId(solution.forces.targetNaifId);
			propagator.setInactiveInteriorMassPolicy(solution.forces.inactiveInteriorMassPolicy);
			if (auto* ephemeris = getDefaultEphemerisSource())
				propagator = propagator.forEncounterMonitoring(*ephemeris);
			propagator.nonGrav() = solution.nonGrav;
			return propagator;
		}

		std::vector<NBodyPropagator::StateVel> receptionStates(const EphemerisSolution& solution,
		                                                       const std::vector<double>& mjdTdb, bool nBody,
		                                                       const NBodyPropagator& propagator,
		                                                       std::stop_token cancel)
		{
			if (nBody) {
				return propagator.propagateAllWithVelocityMjd(solution.position, solution.velocity,
				                                              solution.epochMjdTdb, mjdTdb, solution.nonGrav, cancel);
			}

			std::vector<NBodyPropagator::StateVel> states(mjdTdb.size());
			const StateVectors initial(solution.position, solution.velocity / Constants::GAUSS_K);
			for (std::size_t i = 0; i < mjdTdb.size(); ++i) {
				if (cancel.stop_requested())
					return {};
				const StateVectors state =
				    Calculations::propagateKeplerState(initial, mjdTdb[i] - solution.epochMjdTdb);
				if (!state.isValid())
					throw std::runtime_error("Two-body ephemeris propagation failed");
				states[i] = {state.first, state.second * Constants::GAUSS_K};
			}
			return states;
		}

		Vec3D toPlane(const Vec3D& ecliptic, EphemerisPlane plane, double cosObliquity, double sinObliquity) noexcept
		{
			if (plane == EphemerisPlane::Ecliptic)
				return ecliptic;
			return Vec3D(ecliptic.x, ecliptic.y * cosObliquity - ecliptic.z * sinObliquity,
			             ecliptic.y * sinObliquity + ecliptic.z * cosObliquity);
		}

		double angleBetween(const Vec3D& a, const Vec3D& b) noexcept
		{
			const double norms = a.magnitude() * b.magnitude();
			if (!(norms > 0.0))
				return std::numeric_limits<double>::quiet_NaN();
			return std::acos(std::clamp(a.dot(b) / norms, -1.0, 1.0));
		}

		Vec3D directionInMovingFrame(const Vec3D& direction, const Vec3D& velocity)
		{
			const double length = direction.magnitude();
			if (!(length > 0.0))
				return direction;
			const Vec3D n = direction / length;
			const Vec3D beta = velocity / Constants::AU_PER_DAY;
			const double beta2 = beta.dot(beta);
			if (!(beta2 < 1.0))
				return direction;
			const double gamma = 1.0 / std::sqrt(1.0 - beta2);
			const double dot = n.dot(beta);
			return (n / gamma + beta * (gamma * dot / (gamma + 1.0) - 1.0)) / (1.0 - dot);
		}

		std::optional<Eigen::MatrixXd> covarianceFactor(const EphemerisSolution& solution)
		{
			const auto& covariance = solution.covariance;
			const int dimension = 6 + solution.nonGrav.fitParameterCount();
			if (covariance.rows() != dimension || covariance.cols() != dimension || !covariance.allFinite() ||
			    (covariance.diagonal().array() <= 0.0).any())
				return std::nullopt;
			const Eigen::VectorXd scale = covariance.diagonal().cwiseSqrt();
			const Eigen::MatrixXd correlation =
			    scale.cwiseInverse().asDiagonal() * covariance * scale.cwiseInverse().asDiagonal();
			if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-10))
				return std::nullopt;
			const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(correlation);
			if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() < -1e-10)
				return std::nullopt;
			return scale.asDiagonal() * eigen.eigenvectors() * eigen.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
		}

		std::optional<EphemerisUncertainty> skyUncertainty(const Vec3D& position, const Eigen::MatrixXd& weighted,
		                                                    EphemerisPlane plane, double c, double s)
		{
			const double radius = position.magnitude();
			const double horizontal = std::hypot(position.x, position.y);
			if (!(radius > 0.0) || !(horizontal > 0.0) || !weighted.allFinite())
				return std::nullopt;
			Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
			if (plane == EphemerisPlane::Equatorial)
				rotation << 1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c;
			const Eigen::MatrixXd projected = rotation * weighted;
			const Eigen::RowVector3d east(-position.y / horizontal, position.x / horizontal, 0.0);
			const Eigen::RowVector3d north(-position.x * position.z / (radius * horizontal),
			                               -position.y * position.z / (radius * horizontal), horizontal / radius);
			Eigen::MatrixXd sky(2, weighted.cols());
			sky.row(0) = east * projected / radius;
			sky.row(1) = north * projected / radius;
			const Eigen::Matrix2d covariance = sky * sky.transpose();
			const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(covariance);
			if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() < -1e-20)
				return std::nullopt;
			constexpr double radiansToArcseconds = 206264.80624709636;
			EphemerisUncertainty out;
			out.longitude3SigmaArcsec = 3.0 * std::sqrt(std::max(0.0, covariance(0, 0))) * radiansToArcseconds;
			out.latitude3SigmaArcsec = 3.0 * std::sqrt(std::max(0.0, covariance(1, 1))) * radiansToArcseconds;
			out.major3SigmaArcsec = 3.0 * std::sqrt(std::max(0.0, eigen.eigenvalues()[1])) * radiansToArcseconds;
			out.minor3SigmaArcsec = 3.0 * std::sqrt(std::max(0.0, eigen.eigenvalues()[0])) * radiansToArcseconds;
			const auto major = eigen.eigenvectors().col(1);
			out.positionAngleDeg = (eigen.eigenvalues()[1] - eigen.eigenvalues()[0]) > 1e-12 * eigen.eigenvalues()[1]
			                           ? std::remainder(std::atan2(major[1], major[0]) * Constants::RAD_TO_DEG, 180.0)
			                           : std::numeric_limits<double>::quiet_NaN();
			out.area3SigmaArcsec2 = Constants::PI * out.major3SigmaArcsec * out.minor3SigmaArcsec;
			out.sky3SigmaArcsec = std::hypot(out.major3SigmaArcsec, out.minor3SigmaArcsec);
			const Eigen::RowVector3d radial(position.x / radius, position.y / radius, position.z / radius);
			out.range3SigmaKm = 3.0 * (radial * projected).norm() * (Constants::AU_IN_METERS / 1000.0);
			return out;
		}

	}  // namespace

	std::size_t ephemerisEpochCount(double startMjd, double endMjd, double stepDays)
	{
		if (!std::isfinite(startMjd) || !std::isfinite(endMjd) || !std::isfinite(stepDays))
			throw std::invalid_argument("Ephemeris epochs must be finite");
		if (!(stepDays > 0.0))
			throw std::invalid_argument("Ephemeris step must be positive");
		if (endMjd < startMjd)
			throw std::invalid_argument("Ephemeris end precedes its start");
		if (startMjd == endMjd)
			return 1;
		const double largestEpoch = std::max(std::abs(startMjd), std::abs(endMjd));
		const double resolution = std::nextafter(largestEpoch, std::numeric_limits<double>::infinity()) - largestEpoch;
		if (stepDays < resolution)
			throw std::invalid_argument("Ephemeris step is below epoch resolution");
		std::size_t first = 1, end = kMaxEphemerisEpochs + 1;
		while (first < end) {
			const std::size_t middle = first + (end - first) / 2;
			if (startMjd + static_cast<double>(middle) * stepDays <= endMjd)
				first = middle + 1;
			else
				end = middle;
		}
		return first;
	}

	EphemerisResult computeEphemeris(const EphemerisRequest& request, const FitRunControl& control)
	{
		const auto& solution = request.solution;
		if (!finite(solution.position) || !finite(solution.velocity) || !std::isfinite(solution.epochMjdTdb))
			throw std::invalid_argument("Ephemeris requires a finite solution state and epoch");
		const std::size_t count = ephemerisEpochCount(request.startMjd, request.endMjd, request.stepDays);
		if (count > kMaxEphemerisEpochs)
			throw std::invalid_argument("Too many ephemeris epochs");

		const ObserverSite site = resolveObserver(request);
		const bool nBody = solution.forces.enabled || solution.nonGrav.enabled();
		if (nBody) {
			auto* ephemeris = getDefaultEphemerisSource();
			if (!ephemeris || !ephemeris->isAvailable())
				throw std::runtime_error("Planetary ephemerides are unavailable");
			if (solution.forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
				throw std::runtime_error("SB441 was requested but is unavailable");
		}

		EphemerisResult result;
		result.nBody = nBody;
		result.nonGravitational = solution.nonGrav.enabled();
		result.lightTimeApplied = request.lightTime && request.center != EphemerisCenter::Heliocentric;
		const auto cancelled = [&] {
			result.cancelled = true;
			return result;
		};
		if (control.cancel.stop_requested())
			return cancelled();

		std::vector<EphemerisRow> rows(count);
		std::vector<double> mjdTdb(count);
		for (std::size_t i = 0; i < count; ++i) {
			if (control.cancel.stop_requested())
				return cancelled();
			auto& row = rows[i];
			row.mjd = request.startMjd + static_cast<double>(i) * request.stepDays;
			row.mjdTdb = time::toMjdTdb(row.mjd, request.scale);
			row.mjdTt = request.scale == time::Scale::TT ? row.mjd : time::fromMjdTdb(row.mjdTdb, time::Scale::TT);
			row.mjdUtc = request.scale == time::Scale::UTC ? row.mjd : time::fromMjdTdb(row.mjdTdb, time::Scale::UTC);
			mjdTdb[i] = row.mjdTdb;
		}

		if (control.cancel.stop_requested())
			return cancelled();
		if (control.progress)
			control.progress("Propagating ephemeris");
		const auto propagator = ephemerisPropagator(solution);
		const auto factor = request.uncertainty ? covarianceFactor(solution) : std::nullopt;
		result.uncertaintyComputed = factor.has_value();
		std::vector<NBodyPropagator::PropagationResultFull> variational;
		std::vector<NBodyPropagator::StateVel> states;
		std::vector<double> stateEpochs = mjdTdb;
		std::vector<double> emissionOffsets(count, 0.0);
		if (result.lightTimeApplied) {
			const auto reception = receptionStates(solution, mjdTdb, nBody, propagator, control.cancel);
			if (control.cancel.stop_requested())
				return cancelled();
			EphemeridesContainer earth(EphemeridesContainer::Planet::Earth);
			std::vector<Vec3D> observers(count);
			for (std::size_t i = 0; i < count; ++i) {
				const auto earthState = earth.getMjd(rows[i].mjdTdb);
				if (earth.error() != EphemeridesContainer::Error::NoError)
					throw std::runtime_error(earth.errorString());
				obs::ObserverState observer;
				if (site.ground)
					observer = obs::groundObserverStateEclAuMjd(rows[i].mjdTdb, rows[i].mjdUtc, site.data);
				observers[i] = earthState.R + observer.offsetEclAu;
			}
			const auto light = obs::solveLightTimeSeries(solution.position, solution.velocity, solution.epochMjdTdb,
			                                              mjdTdb, reception, observers, std::vector<bool>(count, true),
			                                              nBody, propagator, solution.nonGrav, control.cancel);
			if (control.cancel.stop_requested())
				return cancelled();
			stateEpochs = light.emissionMjd;
			emissionOffsets = light.emissionOffsetDays;
			states = light.states;
		}
		if (factor && nBody) {
			variational = propagator.propagateAllWithStmAndSensMjd(solution.position, solution.velocity,
			    solution.epochMjdTdb, stateEpochs, solution.nonGrav, nullptr, control.cancel);
			states.clear();
			states.reserve(variational.size());
			for (std::size_t i = 0; i < variational.size(); ++i)
				states.push_back({variational[i].R + variational[i].V * emissionOffsets[i], variational[i].V});
		} else if (!result.lightTimeApplied) {
			states = receptionStates(solution, mjdTdb, nBody, propagator, control.cancel);
		}
		if (control.cancel.stop_requested())
			return cancelled();
		if (states.size() != count)
			throw std::runtime_error("Incomplete ephemeris propagation");

		EphemeridesContainer earth(EphemeridesContainer::Planet::Earth);
		const double cosObliquity = std::cos(Constants::EARTH_INCL);
		const double sinObliquity = std::sin(Constants::EARTH_INCL);
		for (std::size_t i = 0; i < count; ++i) {
			if (control.cancel.stop_requested())
				return cancelled();
			auto& row = rows[i];
			Vec3D target = states[i].R;
			Vec3D targetVelocity = states[i].V;
			Vec3D heliocentricTarget = target;
			Vec3D heliocentricVelocity = targetVelocity;
			Vec3D relative = target;
			Vec3D relativeVelocity = targetVelocity;
			Vec3D rangeVelocity = targetVelocity;
			Vec3D sunToTarget = heliocentricTarget;
			Eigen::Matrix<double, 3, 6> lightTimeJacobian;
			lightTimeJacobian << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero();
			Eigen::Matrix<double, 3, 6> velocityJacobian;
			velocityJacobian << Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();

			if (request.center != EphemerisCenter::Heliocentric) {
				const auto earthState = earth.getMjd(row.mjdTdb);
				if (earth.error() != EphemeridesContainer::Error::NoError)
					throw std::runtime_error(earth.errorString());
				obs::ObserverState observer;
				if (site.ground)
					observer = obs::groundObserverStateEclAuMjd(row.mjdTdb, row.mjdUtc, site.data);
				const Vec3D observerVelocity = earthState.V + observer.velocityEclAuPerDay;
				rangeVelocity = targetVelocity - observerVelocity;

				if (result.lightTimeApplied) {
					const Vec3D observerPosition = earthState.R + observer.offsetEclAu;
					const Vec3D emissionVelocity = targetVelocity;
					const Vec3D emissionSunVelocity = obs::sunVelocityMjd(stateEpochs[i]);
					const Vec3D receptionSunVelocity = obs::sunVelocityMjd(row.mjdTdb);
					target += obs::sunPositionMjd(stateEpochs[i]) + emissionSunVelocity * emissionOffsets[i] -
					          obs::sunPositionMjd(row.mjdTdb);
					const Vec3D ray = target - observerPosition;
					const Vec3D direction = ray / ray.magnitude();
					const Vec3D targetBarycentricVelocity = emissionVelocity + emissionSunVelocity;
					const Vec3D observerBarycentricVelocity = observerVelocity + receptionSunVelocity;
					const double emissionRate =
					    (Constants::AU_PER_DAY + direction.dot(observerBarycentricVelocity)) /
					    (Constants::AU_PER_DAY + direction.dot(targetBarycentricVelocity));
					rangeVelocity = targetBarycentricVelocity - observerBarycentricVelocity;
					targetVelocity = targetBarycentricVelocity * emissionRate - receptionSunVelocity;
					if (factor) {
						lightTimeJacobian.leftCols<3>() = obs::emissionPositionFactor(ray, targetBarycentricVelocity);
						const Vec3D acceleration =
						    (nBody ? propagator.accelerationAtMjd(heliocentricTarget, emissionVelocity, stateEpochs[i], solution.nonGrav)
						           : heliocentricTarget * (-Constants::GM_SUN / std::pow(heliocentricTarget.magnitude(), 3))) +
						    obs::sunAccelerationMjd(stateEpochs[i]);
						const Eigen::Vector3d a(acceleration.x, acceleration.y, acceleration.z);
						const Eigen::RowVector3d n(direction.x, direction.y, direction.z);
						velocityJacobian.leftCols<3>() = -a * n / (Constants::AU_PER_DAY + direction.dot(targetBarycentricVelocity));
					}
					// Trace sunlight to the target at its emission epoch, then account for
					// aberration of that direction in the target's barycentric rest frame.
					const double heliocentricDistance = heliocentricTarget.magnitude();
					if (!(heliocentricDistance > 0.0))
						throw std::runtime_error("Solar distance is zero");
					double sunDelay = heliocentricDistance / Constants::AU_PER_DAY;
					for (int pass = 0; pass < 2; ++pass) {
						sunToTarget =
						    heliocentricTarget - obs::solarDisplacementDuringDelay(time::jdFromMjd(stateEpochs[i]), sunDelay);
						sunDelay = sunToTarget.magnitude() / Constants::AU_PER_DAY;
					}
					row.heliocentricRateAuPerDay = heliocentricTarget.dot(emissionVelocity) / heliocentricDistance;
					row.phase = angleBetween(directionInMovingFrame(sunToTarget, targetBarycentricVelocity),
					                         target - observerPosition);
				}

				relative = target - earthState.R - observer.offsetEclAu;
				relativeVelocity = targetVelocity - earthState.V - observer.velocityEclAuPerDay;
				const Vec3D observerHelio = earthState.R + observer.offsetEclAu;
				if (!result.lightTimeApplied)
					row.phase = angleBetween(heliocentricTarget, relative);
				if (result.lightTimeApplied) {
					Vec3D sunVelocity;
					const double jd = time::jdFromMjd(row.mjdTdb);
					const Vec3D solarRay = observerHelio * -1.0 +
					                       obs::solarDisplacementDuringDelay(jd, observerHelio.magnitude() /
					                                                                Constants::AU_PER_DAY);
					if (auto* ephemeris = getDefaultEphemerisSource())
						if (const auto sun = ephemeris->barycentricSunMjd(row.mjdTdb))
							sunVelocity = sun->V;
					const Vec3D observerBarycentricVelocity = observerVelocity + sunVelocity;
					const Vec3D deflected = obs::absoluteDeflectedDirection(relative, heliocentricTarget, observerHelio);
					row.elongation = angleBetween(directionInMovingFrame(solarRay, observerBarycentricVelocity * -1.0),
					                              directionInMovingFrame(deflected, observerBarycentricVelocity * -1.0));
				} else {
					row.elongation = angleBetween(observerHelio * -1.0, relative);
				}
			}

			row.position = toPlane(relative, request.plane, cosObliquity, sinObliquity);
			row.velocity = toPlane(relativeVelocity, request.plane, cosObliquity, sinObliquity);
			if (factor) {
				Eigen::MatrixXd transport(6, factor->rows());
				if (nBody) {
					transport.leftCols<6>() = variational[i].Phi;
					for (int p = 0; p < solution.nonGrav.fitParameterCount(); ++p)
						transport.col(6 + p) = variational[i].SensA.col(solution.nonGrav.sensitivityColumnForFitParameter(p));
				} else {
					const StateVectors initial(solution.position, solution.velocity / Constants::GAUSS_K);
					for (int p = 0; p < 6; ++p) {
						const Vec3D unit(p % 3 == 0 ? 1.0 : 0.0, p % 3 == 1 ? 1.0 : 0.0, p % 3 == 2 ? 1.0 : 0.0);
						const StateVectors tangent(p < 3 ? unit : Vec3D{}, p >= 3 ? unit / Constants::GAUSS_K : Vec3D{});
						const auto derivative = Calculations::propagateKeplerStateRate(
						    initial, stateEpochs[i] - solution.epochMjdTdb, tangent, 0.0);
						transport.col(p) << derivative.first.x, derivative.first.y, derivative.first.z,
						    derivative.second.x * Constants::GAUSS_K, derivative.second.y * Constants::GAUSS_K,
						    derivative.second.z * Constants::GAUSS_K;
					}
				}
				row.uncertainty = skyUncertainty(row.position, lightTimeJacobian * transport * *factor,
				                                 request.plane, cosObliquity, sinObliquity);
				if (row.uncertainty) {
					const double distance = relative.magnitude();
					const Eigen::Vector3d direction(relative.x / distance, relative.y / distance, relative.z / distance);
					const Eigen::Vector3d velocity(rangeVelocity.x, rangeVelocity.y, rangeVelocity.z);
					// Differentiate n dot v, including both the line of sight and emission velocity.
					const Eigen::RowVector3d positionGradient =
					    (velocity - direction * direction.dot(velocity)).transpose() / distance;
					const Eigen::RowVectorXd gradient =
					    (positionGradient * lightTimeJacobian + direction.transpose() * velocityJacobian) * transport * *factor;
					row.uncertainty->rangeRate3SigmaKmPerSecond =
					    3.0 * gradient.norm() * (Constants::AU_IN_METERS / 1000.0) / 86400.0;
				}
			}
			row.distanceAu = row.position.magnitude();
			row.heliocentricAu = sunToTarget.magnitude();
			row.distanceRateAuPerDay = row.distanceAu > 0.0 ? relative.dot(rangeVelocity) / row.distanceAu : 0.0;
			if (!result.lightTimeApplied)
				row.heliocentricRateAuPerDay =
				    row.heliocentricAu > 0.0 ? heliocentricTarget.dot(heliocentricVelocity) / row.heliocentricAu : 0.0;
			row.oneWayLightTimeDays = row.distanceAu / Constants::AU_PER_DAY;

			const Vec3D& p = row.position;
			const Vec3D& v = row.velocity;
			double longitude = std::atan2(p.y, p.x);
			if (longitude < 0.0)
				longitude += 2.0 * Constants::PI;
			row.longitude = longitude;
			row.latitude = std::asin(std::clamp(p.z / row.distanceAu, -1.0, 1.0));
			const double rho2 = p.x * p.x + p.y * p.y;
			if (rho2 > 0.0) {
				const double rho = std::sqrt(rho2);
				const double r2 = row.distanceAu * row.distanceAu;
				row.longitudeRate = (p.x * v.y - p.y * v.x) / rho2 * std::cos(row.latitude);
				row.latitudeRate = (v.z * r2 - p.z * p.dot(v)) / (r2 * rho);
			}

			if (solution.absoluteMagnitude && request.center != EphemerisCenter::Heliocentric &&
			    solution.magnitudeType == EphemerisMagnitudeType::Asteroid)
				row.visualMagnitude = photometry::visualMagnitude(*solution.absoluteMagnitude, solution.slopeG,
				                                                  row.heliocentricAu, row.distanceAu, row.phase);
			else if (solution.absoluteMagnitude && request.center != EphemerisCenter::Heliocentric)
				row.visualMagnitude =
				    photometry::cometMagnitude(*solution.absoluteMagnitude, solution.cometSlope, row.heliocentricAu,
				                               row.distanceAu, row.phase, solution.cometPhaseCoefficient);
		}

		if (control.cancel.stop_requested())
			return cancelled();
		result.rows = std::move(rows);
		return result;
	}

}  // namespace orbital_core
