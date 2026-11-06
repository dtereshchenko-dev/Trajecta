#include "hazard/CloseApproachTargetPlaneDifferential.h"
#include <span>
#include <utility>

#include "hazard/CloseApproachKinematics.h"
#include "NBodyPropagator.h"
#include "hazard/EncounterPlane.h"
#include "hazard/HazardAudit.h"
#include "orbital_core/Constants.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/nongrav/NonGravProfile.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace orbital_core::hazard {
	namespace {
		using VecX = Eigen::VectorXd;
		using MatX = Eigen::MatrixXd;

		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;

		[[nodiscard]] Vec3D targetAccelerationAtMjd(IEphemerisSource& eph, Planet planet, double mjd)
		{
			constexpr double intervalDays = 0.001;
			const auto id = toEphemerisPlanet(planet);
			return (eph.queryMjd(id, mjd + intervalDays).V - eph.queryMjd(id, mjd - intervalDays).V) /
			       (2.0 * intervalDays);
		}

		struct PropagatedFrame {
			MatX dStateDParams;
			Vec3D rTarget;
			Vec3D vTarget;
			bool ok{false};
		};

		[[nodiscard]] bool isInteriorEncounter(const CloseApproachEvent& event) noexcept
		{
			const double radiusKm = planetRadiusKm(event.planet);
			return radiusKm > 0.0 && event.distAu >= 0.0 && event.distAu * kAuKm < radiusKm;
		}

		[[nodiscard]] NonGravModel nonGravFromParameters(const NonGravModel& nominal, const VecX& params)
		{
			NonGravModel ng = nominal;
			if (ng.numFitParams() == 0)
				return ng;

			std::vector<double> p(static_cast<std::size_t>(params.size()));
			for (Eigen::Index i = 0; i < params.size(); ++i)
				p[static_cast<std::size_t>(i)] = params(i);

			NonGravProfile::updateFromParams(ng, p, 6);
			return ng;
		}

		[[nodiscard]] bool validDifferentialInput(const CloseApproachTargetPlaneDifferentialInput& input)
		{
			return input.propagator && input.nominalNonGrav && input.eph && input.parameters.size() >= 6 &&
			       input.parameterCovariance.rows() == input.parameters.size() &&
			       input.parameterCovariance.cols() == input.parameters.size() &&
			       input.parameterStepPerSigma.size() == input.parameters.size();
		}

		[[nodiscard]] std::vector<PropagatedFrame> propagateFrames(
		    const CloseApproachTargetPlaneDifferentialInput& input, const NonGravModel& ng,
		    std::span<const double> targetJds)
		{
			std::vector<PropagatedFrame> out(targetJds.size());
			if (targetJds.empty())
				return out;

			const int nParams = static_cast<int>(input.parameters.size());
			for (PropagatedFrame& frame : out)
				frame.dStateDParams = MatX::Zero(6, nParams);

			try {
				const Vec3D r0(input.parameters(0), input.parameters(1), input.parameters(2));
				const Vec3D v0(input.parameters(3), input.parameters(4), input.parameters(5));
				const std::vector<double> times(targetJds.begin(), targetJds.end());
				const bool continuousEncounter = input.propagator->encounterStepPointCount() > 0;

				if (nParams > 6 && ng.numFitParams() > 0) {
					const auto states = continuousEncounter
					                        ? input.propagator->propagateAllWithStmAndSens(
					                              r0, v0, input.epochTdb, times, ng, input.integrationReport)
					                        : input.propagator->propagateAllWithStmAndSensChained(
					                              r0, v0, input.epochTdb, times, ng, input.integrationReport);
					if (states.size() != targetJds.size())
						return out;

					for (std::size_t i = 0; i < states.size(); ++i) {
						PropagatedFrame& frame = out[i];
						frame.rTarget = states[i].R;
						frame.vTarget = states[i].V;
						frame.dStateDParams.block(0, 0, 6, 6) = states[i].Phi;
						if (ng.model() == NgModel::Yarkovsky) {
							frame.dStateDParams.col(6) = states[i].SensA.col(1);
						} else {
							const int nExtra = std::min(nParams - 6, ng.numFitParams());
							frame.dStateDParams.block(0, 6, 6, nExtra) = states[i].SensA.leftCols(nExtra);
						}
						frame.ok = true;
					}
				} else {
					const auto states = continuousEncounter
					                        ? input.propagator->propagateAllWithStm(r0, v0, input.epochTdb, times, ng,
					                                                                input.integrationReport)
					                        : input.propagator->propagateAllWithStmChained(
					                              r0, v0, input.epochTdb, times, ng, input.integrationReport);
					if (states.size() != targetJds.size())
						return out;

					for (std::size_t i = 0; i < states.size(); ++i) {
						PropagatedFrame& frame = out[i];
						frame.rTarget = states[i].R;
						frame.vTarget = states[i].V;
						frame.dStateDParams.block(0, 0, 6, 6) = states[i].Phi;
						frame.ok = true;
					}
				}
			} catch (...) {
				for (PropagatedFrame& frame : out)
					frame.ok = false;
			}
			return out;
		}

		[[nodiscard]] bool alignFrameMinimum(const CloseApproachTargetPlaneDifferentialInput& input,
		                                     CloseApproachEvent& event, PropagatedFrame& frame)
		{
			constexpr double rdotTolerance = 1.0e-13;
			const auto planetId = toEphemerisPlanet(event.planet);
			auto radial = [](const PropagatedFrame& f, const EphemerisRecord& planet) {
				const Vec3D r = f.rTarget - planet.R;
				return r.dot(f.vTarget - planet.V) / r.magnitude();
			};
			const auto planet = input.eph->queryJD(planetId, event.jd);
			const double initialRdot = radial(frame, planet);
			if (std::abs(initialRdot) <= rdotTolerance)
				return true;

			auto segment = input;
			segment.epochTdb = event.jd;
			segment.parameters.head(6) << frame.rTarget.x, frame.rTarget.y, frame.rTarget.z, frame.vTarget.x,
			    frame.vTarget.y, frame.vTarget.z;
			const auto ng = nonGravFromParameters(*input.nominalNonGrav, input.parameters);
			const PropagatedFrame initial = frame;
			auto at = [&](double jd) -> std::optional<PropagatedFrame> {
				auto frames = propagateFrames(segment, ng, std::span<const double>(&jd, 1));
				if (frames.empty() || !frames.front().ok)
					return std::nullopt;
				auto out = std::move(frames.front());
				MatX derivatives = out.dStateDParams.leftCols(6) * initial.dStateDParams;
				if (derivatives.cols() > 6)
					derivatives.rightCols(derivatives.cols() - 6) +=
					    out.dStateDParams.rightCols(derivatives.cols() - 6);
				out.dStateDParams = std::move(derivatives);
				return out;
			};
			const double timeScale = (frame.rTarget - planet.R).magnitude() / (frame.vTarget - planet.V).magnitude();
			if (!(timeScale > 0.0) || !std::isfinite(initialRdot))
				return false;
			double span = std::max(1.0e-6, timeScale / std::max(1, input.propagator->encounterStepPointCount()));
			double a = event.jd, b = event.jd, fa = initialRdot, fb = initialRdot;
			for (int it = 0; it < 12 && !(fa < 0.0 && fb > 0.0); ++it) {
				const double trial = event.jd + (initialRdot > 0.0 ? -span : span);
				const auto f = at(trial);
				if (!f)
					return false;
				const double value = radial(*f, input.eph->queryJD(planetId, trial));
				if (!std::isfinite(value))
					return false;
				if (initialRdot > 0.0) {
					a = trial;
					fa = value;
				} else {
					b = trial;
					fb = value;
				}
				span *= 2.0;
			}
			if (!(fa < 0.0 && fb > 0.0))
				return false;
			double previous = event.jd;
			for (int it = 0; it < 50; ++it) {
				const double jd = b - fb * (b - a) / (fb - fa);
				auto f = at(jd);
				if (!f)
					return false;
				const auto p = input.eph->queryJD(planetId, jd);
				const double value = radial(*f, p);
				if (!std::isfinite(value))
					return false;
				if (std::abs(value) <= rdotTolerance || jd == previous || jd == a || jd == b) {
					audit::write("tp-stm-minimum", [&](auto& log) {
						log << " scanJd=" << event.jd << " jd=" << jd << " initialRdot=" << initialRdot
						    << " finalRdot=" << value << " iterations=" << it + 1;
					});
					frame = std::move(*f);
					event.jd = jd;
					event.rPlanet = p.R;
					event.vPlanet = p.V;
					return true;
				}
				if (value < 0.0) {
					a = jd;
					fa = value;
				} else {
					b = jd;
					fb = value;
				}
				previous = jd;
			}
			return false;
		}

		[[nodiscard]] bool projectFrameMinimum(const CloseApproachTargetPlaneDifferentialInput& input,
		                                       const CloseApproachEvent& event, PropagatedFrame& frame)
		{
			const Vec3D r = frame.rTarget - event.rPlanet;
			const Vec3D v = frame.vTarget - event.vPlanet;
			const double mjd = event.jd - 2400000.5;
			const auto ng = nonGravFromParameters(*input.nominalNonGrav, input.parameters);
			const Vec3D acceleration = input.propagator->accelerationAtMjd(frame.rTarget, frame.vTarget, mjd, ng) -
			                           targetAccelerationAtMjd(*input.eph, event.planet, mjd);
			const double speedSquared = v.dot(v);
			const double radialAcceleration = r.dot(acceleration);
			const double rate = speedSquared + radialAcceleration;
			if (!(rate > 64.0 * std::numeric_limits<double>::epsilon() * (speedSquared + std::abs(radialAcceleration))))
				return false;
			Eigen::Matrix<double, 6, 1> flow;
			flow << v.x, v.y, v.z, acceleration.x, acceleration.y, acceleration.z;
			Eigen::Matrix<double, 1, 6> normal;
			normal << v.x, v.y, v.z, r.x, r.y, r.z;
			const MatX derivatives = frame.dStateDParams - flow * (normal * frame.dStateDParams) / rate;
			if (!derivatives.allFinite())
				return false;
			frame.dStateDParams = derivatives;
			return true;
		}

		[[nodiscard]] std::optional<CloseApproachTargetPlaneTrace> traceFromFrame(
		    const CloseApproachTargetPlaneDifferentialInput& input, CloseApproachEvent sourceEvent,
		    PropagatedFrame frame)
		{
			if (!frame.ok)
				return std::nullopt;
			if (!isInteriorEncounter(sourceEvent)) {
				if (!alignFrameMinimum(input, sourceEvent, frame) || !projectFrameMinimum(input, sourceEvent, frame))
					return std::nullopt;
			}

			const EncounterPlanePoint center =
			    encounterPlaneAtHeliocentricState(frame.rTarget, frame.vTarget, sourceEvent.rPlanet,
			                                      sourceEvent.vPlanet, sourceEvent.planet, input.planeKind);
			if (!center.valid)
				return std::nullopt;

			const EncounterPlaneJacobian map =
			    encounterPlaneJacobianAtHeliocentricState(frame.rTarget, frame.vTarget, sourceEvent.rPlanet,
			                                              sourceEvent.vPlanet, sourceEvent.planet, input.planeKind);
			if (!map.valid)
				return std::nullopt;

			const int nParams = static_cast<int>(input.parameters.size());
			const Vec3D rRel = frame.rTarget - sourceEvent.rPlanet;
			const Vec3D vRel = frame.vTarget - sourceEvent.vPlanet;

			audit::write("tp-differential", [&](auto& log) {
				log << " jd=" << sourceEvent.jd << " rx=" << rRel.x << " ry=" << rRel.y << " rz=" << rRel.z
				    << " vx=" << vRel.x << " vy=" << vRel.y << " vz=" << vRel.z;
				for (int i = 0; i < 6; ++i)
					for (int j = 0; j < nParams; ++j)
						log << " stm" << i << '_' << j << '=' << frame.dStateDParams(i, j);
				for (int i = 0; i < 2; ++i)
					for (int j = 0; j < 6; ++j)
						log << " map" << i << '_' << j << '=' << map.dCoordsDTargetState(i, j);
			});

			TargetPlaneDifferential differential;
			differential.coordsKm = center.coordsKm;
			differential.dCoordsDParams = map.dCoordsDTargetState * frame.dStateDParams;
			differential.parameterCovariance = input.parameterCovariance;
			differential.parameterStepPerSigma = input.parameterStepPerSigma;
			differential.bodyRadiusKm = planetRadiusKm(sourceEvent.planet);
			differential.impactCrossSectionBodyRadii = center.gamma;
			differential.sigma = input.sigma;
			differential.jd = sourceEvent.jd;

			TargetPlaneTracePoint trace = makeTargetPlaneTracePoint(differential);
			if (!trace.valid)
				return std::nullopt;

			CloseApproachEvent event = sourceEvent;
			event.rTarget = frame.rTarget;
			event.vTarget = frame.vTarget;
			refreshCloseApproachKinematics(event);
			event.vInfKmS = center.vInfKmS;
			return CloseApproachTargetPlaneTrace{event,        trace,        std::move(differential.dCoordsDParams),
			                                     center.basis, center.gamma, center.vInfKmS};
		}

		[[nodiscard]] bool mapSurfaceFrameToInteriorPericenter(const CloseApproachTargetPlaneDifferentialInput& input,
		                                                       const CloseApproachEvent& sourceEvent,
		                                                       const NonGravModel& ng, double surfaceMjd,
		                                                       PropagatedFrame& frame, CloseApproachEvent& event)
		{
			try {
				const EphemerisPlanet ephPlanet = toEphemerisPlanet(sourceEvent.planet);
				const double surfaceJd = surfaceMjd + 2400000.5;
				const EphemerisRecord surfacePlanet = input.eph->queryMjd(ephPlanet, surfaceMjd);
				audit::write("interior-entry", [&](auto& log) {
					const Vec3D dr = frame.rTarget - surfacePlanet.R;
					const Vec3D dv = frame.vTarget - surfacePlanet.V;
					log << " jd=" << surfaceJd << " mjd=" << surfaceMjd << " rx=" << dr.x << " ry=" << dr.y
					    << " rz=" << dr.z << " vx=" << dv.x << " vy=" << dv.y << " vz=" << dv.z;
					for (int i = 0; i < 6; ++i)
						for (int j = 0; j < frame.dStateDParams.cols(); ++j)
							log << " stm" << i << '_' << j << '=' << frame.dStateDParams(i, j);
				});
				const InteriorPericenterMapping mapping = mapInteriorEncounterToPericenter(
				    frame.rTarget - surfacePlanet.R, frame.vTarget - surfacePlanet.V, sourceEvent.planet);
				const auto pericenterJacobian = mapping.eventJacobian(sourceEvent.planet);
				if (!pericenterJacobian)
					return false;

				const Vec3D r = frame.rTarget - surfacePlanet.R;
				const Vec3D v = frame.vTarget - surfacePlanet.V;
				const double radius = r.magnitude();
				const Vec3D planetAcceleration = targetAccelerationAtMjd(*input.eph, sourceEvent.planet, surfaceMjd);
				const Vec3D incomingAcceleration =
				    input.propagator->accelerationAtMjd(frame.rTarget, frame.vTarget, surfaceMjd, ng) -
				    planetAcceleration;
				const double mu = planetMuKm3PerS2(sourceEvent.planet) * 86400.0 * 86400.0 / (kAuKm * kAuKm * kAuKm);
				const Vec3D outgoingAcceleration = r * (-mu / (radius * radius * radius));
				Eigen::Matrix<double, 6, 1> before, after;
				before << v.x, v.y, v.z, incomingAcceleration.x, incomingAcceleration.y, incomingAcceleration.z;
				after << v.x, v.y, v.z, outgoingAcceleration.x, outgoingAcceleration.y, outgoingAcceleration.z;
				Eigen::Matrix<double, 1, 6> normal;
				normal << r.x / radius, r.y / radius, r.z / radius, 0.0, 0.0, 0.0;
				const double crossingRate = (normal * before).value();
				if (!(crossingRate < -64.0 * std::numeric_limits<double>::epsilon() * v.magnitude()))
					return false;
				const Eigen::Matrix<double, 6, 6> surfaceJacobian =
				    Eigen::Matrix<double, 6, 6>::Identity() + (after - before) * normal / crossingRate;
				const MatX derivatives = *pericenterJacobian * surfaceJacobian * frame.dStateDParams;
				if (!derivatives.allFinite())
					return false;

				const double pericenterMjd = surfaceMjd + mapping.deltaDays;
				const double pericenterJd = 2400000.5 + pericenterMjd;
				audit::write("interior-map", [&](auto& log) {
					log << " jd=" << pericenterJd << " deltaDays=" << mapping.deltaDays;
					for (int i = 0; i < 6; ++i)
						for (int j = 0; j < 6; ++j)
							log << " stm" << i << '_' << j << '=' << mapping.stateTransition(i, j);
				});
				const EphemerisRecord pericenterPlanet = input.eph->queryMjd(ephPlanet, pericenterMjd);
				frame.rTarget = pericenterPlanet.R + mapping.relativePosition;
				frame.vTarget = pericenterPlanet.V + mapping.relativeVelocity;
				frame.dStateDParams = derivatives;

				event = sourceEvent;
				event.jd = pericenterJd;
				event.rPlanet = pericenterPlanet.R;
				event.vPlanet = pericenterPlanet.V;
				event.rTarget = frame.rTarget;
				event.vTarget = frame.vTarget;
				refreshCloseApproachKinematics(event);
				return true;
			} catch (...) {
				return false;
			}
		}
	}  // namespace

	std::optional<CloseApproachTargetPlaneTrace> makeCloseApproachTargetPlaneTracePoint(
	    const CloseApproachTargetPlaneDifferentialInput& input, std::span<const double> intermediateCheckpointJds)
	{
		if (!input.event || !validDifferentialInput(input))
			return std::nullopt;

		const CloseApproachEvent& sourceEvent = *input.event;
		const NonGravModel ng = nonGravFromParameters(*input.nominalNonGrav, input.parameters);
		double targetJd = sourceEvent.jd;
		double surfaceMjd{};
		const bool regularizeInterior = isInteriorEncounter(sourceEvent);
		if (regularizeInterior && targetJd < input.epochTdb)
			return std::nullopt;
		const bool forward = targetJd >= input.epochTdb;
		std::vector<double> propagationJds;
		propagationJds.reserve(intermediateCheckpointJds.size() + 1);
		for (const double jd : intermediateCheckpointJds) {
			if (input.propagator->encounterStepPointCount() > 0)
				break;
			if (!std::isfinite(jd))
				continue;
			const bool isIntermediate =
			    forward ? jd > input.epochTdb && jd < targetJd : jd < input.epochTdb && jd > targetJd;
			if (isIntermediate)
				propagationJds.push_back(jd);
		}
		std::sort(propagationJds.begin(), propagationJds.end());
		propagationJds.erase(std::unique(propagationJds.begin(), propagationJds.end()), propagationJds.end());
		if (!forward)
			std::reverse(propagationJds.begin(), propagationJds.end());
		propagationJds.push_back(targetJd);

		std::vector<PropagatedFrame> frames;
		if (regularizeInterior) {
			try {
				const auto stopped = input.propagator->propagateWithStmToSurface(
				    Vec3D(input.parameters(0), input.parameters(1), input.parameters(2)),
				    Vec3D(input.parameters(3), input.parameters(4), input.parameters(5)), input.epochTdb,
				    propagationJds, toEphemerisPlanet(sourceEvent.planet), planetRadiusKm(sourceEvent.planet),
				    *input.eph, ng, input.integrationReport);
				if (!stopped.hitSurface)
					return std::nullopt;
				targetJd = stopped.jd;
				surfaceMjd = stopped.mjd;
				PropagatedFrame frame;
				frame.rTarget = stopped.state.R;
				frame.vTarget = stopped.state.V;
				const int nParams = static_cast<int>(input.parameters.size());
				frame.dStateDParams = MatX::Zero(6, nParams);
				frame.dStateDParams.leftCols(6) = stopped.state.Phi;
				if (nParams > 6 && ng.numFitParams() > 0) {
					if (ng.model() == NgModel::Yarkovsky)
						frame.dStateDParams.col(6) = stopped.state.SensA.col(1);
					else {
						const int count = std::min(nParams - 6, ng.numFitParams());
						frame.dStateDParams.block(0, 6, 6, count) = stopped.state.SensA.leftCols(count);
					}
				}
				frame.ok = true;
				frames.push_back(std::move(frame));
				audit::write("surface-stm-event", [&](auto& log) {
					const auto planet = input.eph->queryMjd(toEphemerisPlanet(sourceEvent.planet), stopped.mjd);
					log << " jd=" << targetJd << " mjd=" << stopped.mjd
					    << " radiusKm=" << (stopped.state.R - planet.R).magnitude() * kAuKm
					    << " npoint=" << input.propagator->encounterStepPointCount();
					if (input.integrationReport)
						log << " encounterSteps=" << input.integrationReport->encounterSteps
						    << " encounterRestarts=" << input.integrationReport->encounterRestarts
						    << " unconverged=" << input.integrationReport->unconvergedEncounterSteps;
				});
			} catch (const std::exception& error) {
				audit::write("interior-stm-failed", [&](auto& log) { log << " error=" << error.what(); });
				return std::nullopt;
			}
		} else {
			frames = propagateFrames(input, ng, propagationJds);
		}
		if (frames.empty())
			return std::nullopt;

		if (regularizeInterior) {
			const PropagatedFrame& frame = frames.back();
			audit::write("surface-stm-frame", [&](auto& log) {
				log << " jd=" << targetJd << " sourceJd=" << sourceEvent.jd << " valid=" << frame.ok;
				if (frame.ok) {
					log << " rxAu=" << frame.rTarget.x << " ryAu=" << frame.rTarget.y << " rzAu=" << frame.rTarget.z
					    << " vxAuDay=" << frame.vTarget.x << " vyAuDay=" << frame.vTarget.y
					    << " vzAuDay=" << frame.vTarget.z;
				}
			});
		}

		CloseApproachEvent event = sourceEvent;
		if (regularizeInterior &&
		    !mapSurfaceFrameToInteriorPericenter(input, sourceEvent, ng, surfaceMjd, frames.back(), event)) {
			return std::nullopt;
		}
		return traceFromFrame(input, event, frames.back());
	}

	std::vector<std::optional<CloseApproachTargetPlaneTrace>> makeCloseApproachTargetPlaneTracePoints(
	    const CloseApproachTargetPlaneDifferentialInput& input, std::span<const CloseApproachEvent> events)
	{
		std::vector<std::optional<CloseApproachTargetPlaneTrace>> results(events.size());
		if (events.empty() || !validDifferentialInput(input)) {
			return results;
		}

		if (std::any_of(events.begin(), events.end(), isInteriorEncounter)) {
			std::vector<double> checkpoints;
			checkpoints.reserve(events.size());
			for (std::size_t i = 0; i < events.size(); ++i) {
				CloseApproachTargetPlaneDifferentialInput pointInput = input;
				pointInput.event = &events[i];
				results[i] = makeCloseApproachTargetPlaneTracePoint(pointInput, checkpoints);
				checkpoints.push_back(events[i].jd);
			}
			return results;
		}

		const NonGravModel ng = nonGravFromParameters(*input.nominalNonGrav, input.parameters);
		std::vector<double> eventJds;
		eventJds.reserve(events.size());
		for (const CloseApproachEvent& event : events)
			eventJds.push_back(event.jd);
		const std::vector<PropagatedFrame> frames = propagateFrames(input, ng, eventJds);
		for (std::size_t i = 0; i < events.size(); ++i)
			results[i] = traceFromFrame(input, events[i], frames[i]);
		return results;
	}

}  // namespace orbital_core::hazard
