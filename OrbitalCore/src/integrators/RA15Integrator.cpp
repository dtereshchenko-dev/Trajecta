#include "RA15Integrator.h"
#include "RA15DenseOutput.h"
#include <optional>
#include <span>
#include "orbital_core/Constants.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "../propagator/VariationalEqs.h"
#include "RadauCore.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

	using namespace orbital_core::integrators::radau;

	inline void addCompensated(double& sum, double& compensation, double increment) noexcept
	{
		const double corrected = increment - compensation;
		const double total = sum + corrected;
		compensation = (total - sum) - corrected;
		sum = total;
	}

	constexpr int kPhysicalComponents = 3;
	constexpr int kStmEnd = kPhysicalComponents * (1 + 6);
	constexpr int kAccelerationSensitivityEnd = kStmEnd + kPhysicalComponents * 3;

	constexpr double controlScale(int component)
	{
		return component < kStmEnd                       ? 1.0
		       : component < kAccelerationSensitivityEnd ? Constants::GM_SUN
		                                                 : 1.0 / Constants::GAUSS_K;
	}

	template <int N, typename AccelerationEvaluator, typename OutputWriter, typename StepObserver>
	IntegrationReport integrateRadau15(RadauVector<N>& x, RadauVector<N>& v, double tini,
	                                   std::span<const double> targetMjds, std::size_t outputCount,
	                                   const IntegrationOptions& options, AccelerationEvaluator&& evaluate,
	                                   OutputWriter&& writeOutput, StepObserver&& observeStep,
	                                   bool observeAcceptedSteps)
	{
		constexpr double kCorrectorTolerance = 1.0e-16;
		constexpr int kMaximumCorrectorIterations = 12;
		constexpr int kStagnationCheckIteration = 3;
		constexpr double kStepPower = 1.0 / 9.0;

		if (targetMjds.size() != outputCount)
			throw std::invalid_argument("RA15 target/output size mismatch");
		if (!(options.initialStepDays > 0.0) || !std::isfinite(options.initialStepDays)) {
			throw std::invalid_argument("RA15 requires a positive finite initial step");
		}

		IntegrationReport report;
		report.finalState = {Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
		report.reachedMjd = tini;
		const auto cancelledAt = [&](double offset) {
			report.cancelled = true;
			report.finalState = {Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
			report.reachedMjd = tini + offset;
			return report;
		};
		if (options.cancel.stop_requested())
			return cancelledAt(0.0);
		if (targetMjds.empty())
			return report;

		const auto countedEvaluate = [&](const RadauVector<N>& stateX, const RadauVector<N>& stateV, double offsetDays,
		                                 RadauVector<N>& acceleration) {
			++report.rhsEvaluations;
			evaluate(stateX, stateV, offsetDays, acceleration);
		};

		IntegratedStateFull accumulated;
		const auto globalState = [&](const RadauVector<N>& localX, const RadauVector<N>& localV) {
			auto state = radauState(localX, localV);
			if constexpr (N >= 21) {
				if (report.encounterRestarts > 0) {
					if constexpr (N >= 33)
						state.Sens = (state.Phi * accumulated.Sens + state.Sens).eval();
					state.Phi = (state.Phi * accumulated.Phi).eval();
				}
			}
			return state;
		};
		const auto writeGlobalOutput = [&](std::size_t index, double mjd, const RadauVector<N>& localX,
		                                   const RadauVector<N>& localV) {
			if constexpr (N >= 21) {
				if (report.encounterRestarts > 0) {
					auto outputX = localX;
					auto outputV = localV;
					setRadauVariations(outputX, outputV, globalState(localX, localV));
					writeOutput(index, mjd, outputX, outputV);
					return;
				}
			}
			writeOutput(index, mjd, localX, localV);
		};

		double direction = 0.0;
		for (double target : targetMjds) {
			const double delta = target - tini;
			if (std::abs(delta) > 1.0e-13) {
				direction = delta > 0.0 ? 1.0 : -1.0;
				break;
			}
		}

		std::size_t targetIndex = 0;
		while (targetIndex < targetMjds.size() && std::abs(targetMjds[targetIndex] - tini) <= 1.0e-13) {
			writeOutput(targetIndex, targetMjds[targetIndex], x, v);
			++targetIndex;
		}
		report.completedOutputs = targetIndex;
		if (direction == 0.0)
			return report;

		double previousTarget = tini;
		for (double target : targetMjds) {
			if (direction * (target - tini) < -1.0e-12 || direction * (target - previousTarget) < -1.0e-12) {
				throw std::invalid_argument("RA15 targets are not monotonic");
			}
			previousTarget = target;
		}

		const double span = targetMjds.back() - tini;
		constexpr double kAccuracy = RA15Integrator::kAccuracy;

		const RadauConstants& rc = radauConstants();
		double tp = direction * options.initialStepDays;
		if (tp / span > 1.0)
			tp = span;
		else if (tp / span > 0.5)
			tp = 0.5 * span;

		std::array<RadauVector<N>, 8> b;
		std::array<RadauVector<N>, 8> g;
		std::array<RadauVector<N>, 8> e;
		std::array<RadauVector<N>, 8> bd;
		RadauVector<N> f1;
		RadauVector<N> fj;
		RadauVector<N> y;
		RadauVector<N> z;

		int firstStepReductions = 0;
		double lastHv = 0.0;
		auto encounterControl = options.nextEncounterStep;
		auto encounterControlOffset = options.nextEncounterStepOffset;
		std::optional<double> fixedStep;
		double segmentStartOffset = 0.0;
		double timeCompensation = 0.0;
		RadauVector<N> positionCompensation = RadauVector<N>::Zero();
		RadauVector<N> velocityCompensation = RadauVector<N>::Zero();

		for (;;) {
			int completedSteps = 0;
			int failedConvergenceReductions = 0;
			double tm = segmentStartOffset;
			if (options.cancel.stop_requested())
				return cancelledAt(tm);
			bool started = false;
			countedEvaluate(x, v, tm, f1);
			clearRadauStorage(b, e, bd);

			bool restartFirstStep = false;
			bool restartEncounter = false;
			for (;;) {
				if (options.cancel.stop_requested())
					return cancelledAt(tm);
				computeRadauBeta(b, g, rc);
				const double step = tp;
				const double currentMjd = tini + tm;
				if (!std::isfinite(step))
					throw std::runtime_error("RA15 produced a non-finite step");
				if (direction * step <= 0.0)
					throw std::runtime_error("RA15 produced a zero or reversed step");
				if (currentMjd + step == currentMjd)
					throw std::runtime_error("RA15 step is below Julian-date resolution");
				const double step2 = step * step;
				const double absStep = std::abs(step);
				double previousCorrectorError = std::numeric_limits<double>::infinity();
				bool converged = false;

				for (int iteration = 1; iteration <= kMaximumCorrectorIterations; ++iteration) {
					++report.correctorIterations;
					double maximumAcceleration = 0.0;
					double maximumCorrection = 0.0;
					for (int j = 2; j <= 8; ++j) {
						if (options.cancel.stop_requested())
							return cancelledAt(tm);
						const double s = rc.h[j];
						for (int k = 0; k < N; ++k) {
							double temp =
							    rc.w[3] * b[3](k) +
							    s * (rc.w[4] * b[4](k) +
							         s * (rc.w[5] * b[5](k) + s * (rc.w[6] * b[6](k) + s * rc.w[7] * b[7](k))));
							y(k) = x(k) +
							       (s * (step * v(k) + step2 * s *
							                               (f1(k) * rc.w1 + s * (rc.w[1] * b[1](k) +
							                                                     s * (rc.w[2] * b[2](k) + s * temp)))) -
							        positionCompensation(k) - s * step * velocityCompensation(k));

							temp = rc.u[3] * b[3](k) +
							       s * (rc.u[4] * b[4](k) +
							            s * (rc.u[5] * b[5](k) + s * (rc.u[6] * b[6](k) + s * rc.u[7] * b[7](k))));
							z(k) = v(k) +
							       (s * step * (f1(k) + s * (rc.u[1] * b[1](k) + s * (rc.u[2] * b[2](k) + s * temp))) -
							        velocityCompensation(k));
						}

						countedEvaluate(y, z, tm + (s * step - timeCompensation), fj);
						for (int k = 0; k < N; ++k) {
							const double oldG = g[j - 1](k);
							const double gk = (fj(k) - f1(k)) / s;
							switch (j) {
								case 2:
									g[1](k) = gk;
									break;
								case 3:
									g[2](k) = (gk - g[1](k)) * rc.r[1];
									break;
								case 4:
									g[3](k) = ((gk - g[1](k)) * rc.r[2] - g[2](k)) * rc.r[3];
									break;
								case 5:
									g[4](k) = (((gk - g[1](k)) * rc.r[4] - g[2](k)) * rc.r[5] - g[3](k)) * rc.r[6];
									break;
								case 6:
									g[5](k) = ((((gk - g[1](k)) * rc.r[7] - g[2](k)) * rc.r[8] - g[3](k)) * rc.r[9] -
									           g[4](k)) *
									          rc.r[10];
									break;
								case 7:
									g[6](k) =
									    (((((gk - g[1](k)) * rc.r[11] - g[2](k)) * rc.r[12] - g[3](k)) * rc.r[13] -
									      g[4](k)) *
									         rc.r[14] -
									     g[5](k)) *
									    rc.r[15];
									break;
								case 8:
									g[7](k) =
									    ((((((gk - g[1](k)) * rc.r[16] - g[2](k)) * rc.r[17] - g[3](k)) * rc.r[18] -
									       g[4](k)) *
									          rc.r[19] -
									      g[5](k)) *
									         rc.r[20] -
									     g[6](k)) *
									    rc.r[21];
									break;
								default:
									break;
							}

							const double delta = g[j - 1](k) - oldG;
							if (j == 8) {
								const double correction = controlScale(k) * std::abs(delta);
								const double acceleration = controlScale(k) * std::abs(fj(k));
								if (!std::isfinite(correction) || !std::isfinite(acceleration))
									throw std::runtime_error("RA15 corrector produced non-finite coefficients");
								maximumCorrection = std::max(maximumCorrection, correction);
								maximumAcceleration = std::max(maximumAcceleration, acceleration);
							}
							b[j - 1](k) += delta;
							const int firstC = (j - 3) * (j - 2) / 2 + 1;
							for (int l = 1; l <= j - 2; ++l)
								b[l](k) += rc.c[firstC + l - 1] * delta;
						}
					}

					const double correctorError =
					    maximumAcceleration > 0.0
					        ? maximumCorrection / maximumAcceleration
					        : (maximumCorrection == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
					if (std::isfinite(correctorError) &&
					    (correctorError < kCorrectorTolerance ||
					     (iteration >= kStagnationCheckIteration && correctorError >= previousCorrectorError))) {
						double maxB7 = 0.0;
						for (int k = 0; k < kPhysicalComponents; ++k) {
							const double magnitude = controlScale(k) * std::abs(b[7](k));
							if (magnitude > maxB7)
								maxB7 = magnitude;
						}
						if (!fixedStep)
							lastHv = maxB7 * rc.w[7] / std::pow(absStep, 7.0);
						converged = true;
						break;
					}
					previousCorrectorError = correctorError;
				}

				if (!converged && fixedStep)
					++report.unconvergedEncounterSteps;
				if (!converged && !fixedStep) {
					tp *= 0.8;
					clearRadauStorage(b, e, bd);
					++report.rejectedSteps;
					if (++failedConvergenceReductions > 1000)
						throw std::runtime_error("RA15 failed to converge");
					continue;
				}

				if (!started) {
					double next = direction * std::pow(kAccuracy / lastHv, kStepPower);
					if (next / step > 1.4)
						next = step * 1.4;
					if (next / step < 1.0) {
						tp = 0.8 * next;
						++report.rejectedSteps;
						if (++firstStepReductions > 20)
							throw std::runtime_error("RA15 initial step reduction failed");
						restartFirstStep = true;
						break;
					}
					tp = next;
					started = true;
				}

				const RadauVector<N> stepStartX = x;
				const RadauVector<N> stepStartV = v;

				while (targetIndex < targetMjds.size()) {
					const double outputOffset = targetMjds[targetIndex] - tini;
					const double outputDelta = (outputOffset - tm) + timeCompensation;
					if (direction * (outputDelta - step) > 1.0e-10)
						break;
					if (direction * outputDelta < -1.0e-10) {
						throw std::runtime_error("RA15 skipped an output epoch");
					}

					const double s = std::clamp(outputDelta / step, 0.0, 1.0);
					RadauVector<N> outputX;
					RadauVector<N> outputV;
					for (int k = 0; k < N; ++k) {
						double temp = rc.w[3] * b[3](k) +
						              s * (rc.w[4] * b[4](k) +
						                   s * (rc.w[5] * b[5](k) + s * (rc.w[6] * b[6](k) + s * rc.w[7] * b[7](k))));
						outputX(k) =
						    stepStartX(k) +
						    (s * (step * stepStartV(k) +
						          step2 * s *
						              (f1(k) * rc.w1 + s * (rc.w[1] * b[1](k) + s * (rc.w[2] * b[2](k) + s * temp)))) -
						     positionCompensation(k) - s * step * velocityCompensation(k));

						temp = rc.u[3] * b[3](k) +
						       s * (rc.u[4] * b[4](k) +
						            s * (rc.u[5] * b[5](k) + s * (rc.u[6] * b[6](k) + s * rc.u[7] * b[7](k))));
						outputV(k) =
						    stepStartV(k) +
						    (s * step * (f1(k) + s * (rc.u[1] * b[1](k) + s * (rc.u[2] * b[2](k) + s * temp))) -
						     velocityCompensation(k));
					}
					writeGlobalOutput(targetIndex, targetMjds[targetIndex], outputX, outputV);
					++targetIndex;
					report.completedOutputs = targetIndex;
				}

				if constexpr (N == 3) {
					if (observeAcceptedSteps)
						observeStep.prepare(tini, tm, timeCompensation, step, stepStartX, stepStartV,
						                    positionCompensation, velocityCompensation, f1, b);
				}
				for (int k = 0; k < N; ++k) {
					const double positionIncrement =
					    step2 * (f1(k) * rc.w1 + b[1](k) * rc.w[1] + b[2](k) * rc.w[2] + b[3](k) * rc.w[3] +
					             b[4](k) * rc.w[4] + b[5](k) * rc.w[5] + b[6](k) * rc.w[6] + b[7](k) * rc.w[7]);
					const double velocityIncrement =
					    step * (f1(k) + b[1](k) * rc.u[1] + b[2](k) * rc.u[2] + b[3](k) * rc.u[3] + b[4](k) * rc.u[4] +
					            b[5](k) * rc.u[5] + b[6](k) * rc.u[6] + b[7](k) * rc.u[7]);
					addCompensated(x(k), positionCompensation(k), v(k) * step);
					addCompensated(x(k), positionCompensation(k), positionIncrement - velocityCompensation(k) * step);
					addCompensated(v(k), velocityCompensation(k), velocityIncrement);
				}

				addCompensated(tm, timeCompensation, step);
				++completedSteps;
				++report.acceptedSteps;
				if (fixedStep)
					++report.encounterSteps;

				bool continueIntegration = true;
				if constexpr (N >= 21) {
					if (options.variationalObserverOffset) {
						continueIntegration = options.variationalObserverOffset(tini, tm, globalState(x, v));
					} else if (options.variationalObserver) {
						continueIntegration = options.variationalObserver(tini + tm, globalState(x, v));
					}
				}
				if (!continueIntegration ||
				    (observeAcceptedSteps && !observeStep(tini + tm - step, stepStartX, stepStartV, tini + tm, x, v))) {
					report.finalState = {Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
					report.reachedMjd = tini + tm;
					report.stoppedByObserver = true;
					return report;
				}

				if (std::abs((tm - span) - timeCompensation) < 1.0e-10) {
					if (options.cancel.stop_requested())
						return cancelledAt(tm);
					if (targetIndex != targetMjds.size())
						throw std::runtime_error("RA15 did not emit every output");
					report.finalState = {Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
					report.reachedMjd = tini + tm;
					return report;
				}

				countedEvaluate(x, v, tm, f1);

				if (encounterControlOffset || encounterControl) {
					const bool wasEncounter = fixedStep.has_value();
					const IntegratedState acceptedState{Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
					fixedStep = encounterControlOffset ? encounterControlOffset(step, tini, tm, acceptedState)
					                                   : encounterControl(step, tini + tm, acceptedState);
					if (fixedStep && (!std::isfinite(*fixedStep) || direction * *fixedStep <= 0.0))
						throw std::runtime_error("Invalid encounter step");
					if (options.restartOnEncounterExit && wasEncounter && !fixedStep) {
						if constexpr (N >= 21) {
							accumulated = globalState(x, v);
							setRadauVariations(x, v, IntegratedStateFull{});
							for (int k = kPhysicalComponents; k < N; ++k) {
								positionCompensation(k) = 0.0;
								velocityCompensation(k) = 0.0;
							}
						}
						++report.encounterRestarts;
						segmentStartOffset = tm;
						firstStepReductions = 0;
						lastHv = 0.0;
						const double remaining = (span - tm) + timeCompensation;
						tp = direction * options.initialStepDays;
						if (tp / remaining > 1.0)
							tp = remaining;
						else if (tp / remaining > 0.5)
							tp = 0.5 * remaining;
						restartEncounter = true;
						break;
					}
				}
				if (fixedStep) {
					tp = *fixedStep;
				} else {
					tp = direction * std::pow(kAccuracy / lastHv, kStepPower);
					if (tp / step > 1.4)
						tp = step * 1.4;
				}
				const double remaining = (span - tm) + timeCompensation;
				if (direction * (remaining - tp) < 1.0e-8) {
					tp = remaining;
				}

				extrapolateRadauCoefficients(tp / step, b, e, bd, completedSteps);
			}

			if (!restartFirstStep && !restartEncounter)
				throw std::runtime_error("RA15 stopped before reaching target");
		}
	}

}  // namespace

IntegrationReport RA15Integrator::integrateState(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
                                                 const IntegratedState& initial, double epochMjd,
                                                 std::span<const double> targetMjds, std::span<IntegratedState> outputs,
                                                 const IntegrationOptions& options,
                                                 const AcceptedStateStepObserver& observer) const
{
	RadauVector<3> x, v;
	packState(initial.R, initial.V, x, v);

	return integrateRadau15<3>(
	    x, v, epochMjd, targetMjds, outputs.size(), options,
	    [&dynamics, &ng, epochMjd](const RadauVector<3>& stageX, const RadauVector<3>& stageV, double offsetDays,
	                               RadauVector<3>& acceleration) {
		    const Vec3D a = dynamics.accelerationAtOffsetMjd(positionFromRadau(stageX), velocityFromRadau(stageV),
		                                                     epochMjd, offsetDays, ng);
		    acceleration(0) = a.x;
		    acceleration(1) = a.y;
		    acceleration(2) = a.z;
	    },
	    [&outputs](std::size_t index, double, const RadauVector<3>& outputX, const RadauVector<3>& outputV) {
		    outputs[index] = {positionFromRadau(outputX), velocityFromRadau(outputV)};
	    },
	    RA15DenseOutput{observer}, static_cast<bool>(observer));
}

IntegrationReport RA15Integrator::integrateStateWithStm(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
                                                        const IntegratedStateWithStm& initial, double epochMjd,
                                                        std::span<const double> targetMjds,
                                                        std::span<IntegratedStateWithStm> outputs,
                                                        const IntegrationOptions& options) const
{
	RadauVector<21> x, v;
	packState(initial.R, initial.V, x, v);
	for (int column = 0; column < 6; ++column) {
		for (int row = 0; row < 3; ++row) {
			const int index = 3 + 3 * column + row;
			x(index) = initial.Phi(row, column);
			v(index) = initial.Phi(row + 3, column);
		}
	}

	const auto noObserver = [](double, const RadauVector<21>&, const RadauVector<21>&, double, const RadauVector<21>&,
	                           const RadauVector<21>&) { return true; };

	return integrateRadau15<21>(
	    x, v, epochMjd, targetMjds, outputs.size(), options,
	    [&dynamics, &ng, epochMjd](const RadauVector<21>& stageX, const RadauVector<21>& stageV, double offsetDays,
	                               RadauVector<21>& acceleration) {
		    Vec3D a;
		    Eigen::Matrix3d dAdR, dAdV;
		    dynamics.accelerationAndPartialsOffsetMjd(positionFromRadau(stageX), velocityFromRadau(stageV), epochMjd,
		                                              offsetDays, ng, a, dAdR, dAdV);
		    acceleration(0) = a.x;
		    acceleration(1) = a.y;
		    acceleration(2) = a.z;
		    for (int column = 0; column < 6; ++column) {
			    const int base = 3 + 3 * column;
			    for (int row = 0; row < 3; ++row) {
				    double value = 0.0;
				    for (int k = 0; k < 3; ++k) {
					    value += dAdR(row, k) * stageX(base + k) + dAdV(row, k) * stageV(base + k);
				    }
				    acceleration(base + row) = value;
			    }
		    }
	    },
	    [&outputs](std::size_t outputIndex, double, const RadauVector<21>& outputX, const RadauVector<21>& outputV) {
		    auto& output = outputs[outputIndex];
		    output.R = positionFromRadau(outputX);
		    output.V = velocityFromRadau(outputV);
		    for (int column = 0; column < 6; ++column) {
			    for (int row = 0; row < 3; ++row) {
				    const int index = 3 + 3 * column + row;
				    output.Phi(row, column) = outputX(index);
				    output.Phi(row + 3, column) = outputV(index);
			    }
		    }
	    },
	    noObserver, false);
}

IntegrationReport RA15Integrator::integrateStateWithStmAndSens(
    const IIntegrationDynamics& dynamics, const NonGravModel& ng, const IntegratedStateFull& initial, double epochMjd,
    std::span<const double> targetMjds, std::span<IntegratedStateFull> outputs, const IntegrationOptions& options) const
{
	RadauVector<33> x, v;
	packState(initial.R, initial.V, x, v);
	for (int column = 0; column < 6; ++column) {
		for (int row = 0; row < 3; ++row) {
			const int index = 3 + 3 * column + row;
			x(index) = initial.Phi(row, column);
			v(index) = initial.Phi(row + 3, column);
		}
	}
	for (int column = 0; column < 4; ++column) {
		for (int row = 0; row < 3; ++row) {
			const int index = 21 + 3 * column + row;
			x(index) = initial.Sens(row, column);
			v(index) = initial.Sens(row + 3, column);
		}
	}

	const auto noObserver = [](double, const RadauVector<33>&, const RadauVector<33>&, double, const RadauVector<33>&,
	                           const RadauVector<33>&) { return true; };

	return integrateRadau15<33>(
	    x, v, epochMjd, targetMjds, outputs.size(), options,
	    [&dynamics, &ng, epochMjd](const RadauVector<33>& stageX, const RadauVector<33>& stageV, double offsetDays,
	                               RadauVector<33>& acceleration) {
		    const Vec3D r = positionFromRadau(stageX);
		    const Vec3D velocity = velocityFromRadau(stageV);
		    Vec3D a;
		    Eigen::Matrix3d dAdR, dAdV;
		    dynamics.accelerationAndPartialsOffsetMjd(r, velocity, epochMjd, offsetDays, ng, a, dAdR, dAdV);
		    acceleration(0) = a.x;
		    acceleration(1) = a.y;
		    acceleration(2) = a.z;

		    for (int column = 0; column < 6; ++column) {
			    const int base = 3 + 3 * column;
			    for (int row = 0; row < 3; ++row) {
				    double value = 0.0;
				    for (int k = 0; k < 3; ++k) {
					    value += dAdR(row, k) * stageX(base + k) + dAdV(row, k) * stageV(base + k);
				    }
				    acceleration(base + row) = value;
			    }
		    }

		    const Eigen::Matrix<double, 3, 4> forcing =
		        dynamics.nonGravPartialsOffsetMjd(r, velocity, epochMjd, offsetDays, ng);
		    for (int column = 0; column < 4; ++column) {
			    const int base = 21 + 3 * column;
			    for (int row = 0; row < 3; ++row) {
				    double value = forcing(row, column);
				    for (int k = 0; k < 3; ++k) {
					    value += dAdR(row, k) * stageX(base + k) + dAdV(row, k) * stageV(base + k);
				    }
				    acceleration(base + row) = value;
			    }
		    }
	    },
	    [&outputs](std::size_t outputIndex, double, const RadauVector<33>& outputX, const RadauVector<33>& outputV) {
		    auto& output = outputs[outputIndex];
		    output.R = positionFromRadau(outputX);
		    output.V = velocityFromRadau(outputV);
		    for (int column = 0; column < 6; ++column) {
			    for (int row = 0; row < 3; ++row) {
				    const int index = 3 + 3 * column + row;
				    output.Phi(row, column) = outputX(index);
				    output.Phi(row + 3, column) = outputV(index);
			    }
		    }
		    for (int column = 0; column < 4; ++column) {
			    for (int row = 0; row < 3; ++row) {
				    const int index = 21 + 3 * column + row;
				    output.Sens(row, column) = outputX(index);
				    output.Sens(row + 3, column) = outputV(index);
			    }
		    }
	    },
	    noObserver, false);
}
