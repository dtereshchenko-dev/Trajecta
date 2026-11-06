#include "IAS15Integrator.h"
#include <optional>
#include <span>
#include "orbital_core/nongrav/NonGravModel.h"
#include "RadauCore.h"
#include "IAS15Coefficients.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

	using namespace orbital_core::integrators::radau;
	using namespace orbital_core::integrators::ias15;


	constexpr double kSafetyFactor = 0.25;
	constexpr int kPhysicalComponents = 3;

	constexpr double kCorrectorTolerance = 1.0e-16;
	constexpr int kMaxCorrectorIterations = 12;
	constexpr int kStagnationCheckIteration = 3;

	constexpr double kSlowMotionSquared = 1.0e-16;
	constexpr double kMaximumPredictorStepRatio = 20.0;

	[[nodiscard]] inline double seventhRoot(double value) noexcept { return std::pow(value, 1.0 / 7.0); }

	inline void addCompensated(double& sum, double& compensation, double increment) noexcept
	{
		const double y = increment - compensation;
		const double t = sum + y;
		compensation = (t - sum) - y;
		sum = t;
	}

	template <int N>
	[[nodiscard]] double globalRelativeError(const std::array<RadauVector<N>, 8>& b, const RadauVector<N>& acceleration,
	                                         const RadauVector<N>& x, const RadauVector<N>& v, double step) noexcept
	{
		double positionSquared = 0.0;
		double velocitySquared = 0.0;
		for (int k = 0; k < kPhysicalComponents; ++k) {
			positionSquared += x(k) * x(k);
			velocitySquared += v(k) * v(k);
		}
		// REBOUND GLOBAL ignores slowly changing physical accelerations; variations
		// participate in the corrector, but never in this timestep error estimate.
		if (std::abs(velocitySquared * step * step / positionSquared) < kSlowMotionSquared)
			return 0.0;

		double maxAcceleration = 0.0;
		double maxCoefficient = 0.0;
		for (int k = 0; k < kPhysicalComponents; ++k) {
			const double a = std::abs(acceleration(k));
			const double coefficient = std::abs(b[7](k));
			if (std::isnormal(a) && a > maxAcceleration)
				maxAcceleration = a;
			if (std::isnormal(coefficient) && coefficient > maxCoefficient)
				maxCoefficient = coefficient;
		}
		return maxAcceleration > 0.0 ? maxCoefficient / maxAcceleration : 0.0;
	}

	template <int N>
	void computeIas15Beta(const std::array<RadauVector<N>, 8>& b, std::array<RadauVector<N>, 8>& g,
	                      const IAS15Coefficients& rc)
	{
		for (int k = 0; k < N; ++k) {
			for (int row = 1; row <= 7; ++row) {
				double value = 0.0;
				for (int column = 7; column > row; --column)
					value += b[column](k) * rc.d[(column - 2) * (column - 1) / 2 + row];
				g[row](k) = value + b[row](k);
			}
		}
	}

	template <int N>
	void predictIas15State(const RadauVector<N>& x, const RadauVector<N>& v, const RadauVector<N>& csx,
	                       const RadauVector<N>& csv, const std::array<RadauVector<N>, 8>& b, const RadauVector<N>& a0,
	                       double s, double step, RadauVector<N>& y, RadauVector<N>& z)
	{
		for (int k = 0; k < N; ++k) {
			double position = b[7](k);
			double velocity = position;
			for (int degree = 7; degree >= 1; --degree) {
				const double lower = degree == 1 ? a0(k) : b[degree - 1](k);
				position = position * degree * s / (degree + 2) + lower;
				velocity = velocity * degree * s / (degree + 1) + lower;
			}
			y(k) = x(k) + ((position * step * s / 2.0 + v(k)) * step * s - csx(k));
			z(k) = v(k) + (velocity * step * s - csv(k));
		}
	}

	template <int N>
	void predictIas15Coefficients(double ratio, std::array<RadauVector<N>, 8>& b, std::array<RadauVector<N>, 8>& e)
	{
		if (ratio > kMaximumPredictorStepRatio) {
			for (int degree = 1; degree <= 7; ++degree) {
				b[degree].setZero();
				e[degree].setZero();
			}
			return;
		}
		static constexpr auto binomial = [] {
			std::array<std::array<double, 8>, 8> table{};
			for (int n = 0; n <= 7; ++n) {
				table[n][0] = table[n][n] = 1.0;
				for (int k = 1; k < n; ++k)
					table[n][k] = table[n - 1][k - 1] + table[n - 1][k];
			}
			return table;
		}();
		const double q2 = ratio * ratio;
		const double q3 = ratio * q2;
		const double q4 = q2 * q2;
		const std::array<double, 8> powers{0.0, ratio, q2, q3, q4, ratio * q4, q3 * q3, q3 * q4};
		for (int k = 0; k < N; ++k) {
			std::array<double, 8> next{}, correction{};
			for (int degree = 1; degree <= 7; ++degree) {
				correction[degree] = b[degree](k) - e[degree](k);
				for (int source = 7; source >= degree; --source)
					next[degree] += binomial[source][degree] * b[source](k);
			}
			for (int degree = 1; degree <= 7; ++degree) {
				e[degree](k) = powers[degree] * next[degree];
				b[degree](k) = e[degree](k) + correction[degree];
			}
		}
	}

	template <int N, typename AccelerationEvaluator, typename OutputWriter, typename StepObserver>
	IntegrationReport integrateIas15(RadauVector<N>& x, RadauVector<N>& v, double tini,
	                                 std::span<const double> targetMjds, std::size_t outputCount,
	                                 const IntegrationOptions& options, AccelerationEvaluator&& evaluate,
	                                 OutputWriter&& writeOutput, StepObserver&& observeStep, bool observeAcceptedSteps)
	{
		if (targetMjds.size() != outputCount)
			throw std::invalid_argument("IAS15 target/output size mismatch");
		if (!(options.initialStepDays > 0.0) || !std::isfinite(options.initialStepDays)) {
			throw std::invalid_argument("IAS15 requires a positive finite initial step");
		}

		constexpr double epsilon = IAS15Integrator::kEpsilon;

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
				throw std::invalid_argument("IAS15 targets are not monotonic");
			}
			previousTarget = target;
		}

		const double span = targetMjds.back() - tini;
		const auto& rc = kCoefficients;

		double tp = direction * options.initialStepDays;
		if (tp / span > 1.0)
			tp = span;
		else if (tp / span > 0.5)
			tp = 0.5 * span;

		std::array<RadauVector<N>, 8> b;
		std::array<RadauVector<N>, 8> g;
		std::array<RadauVector<N>, 8> e;
		std::array<RadauVector<N>, 8> csb;
		std::array<RadauVector<N>, 8> acceptedB;
		std::array<RadauVector<N>, 8> acceptedE;
		RadauVector<N> f1;
		RadauVector<N> fj;
		RadauVector<N> y;
		RadauVector<N> z;

		RadauVector<N> csx;
		RadauVector<N> csv;

		auto encounterControl = options.nextEncounterStep;
		auto encounterControlOffset = options.nextEncounterStepOffset;
		std::optional<double> fixedStep;
		double segmentStartOffset = 0.0;
		// The physical elapsed time is tm - timeCompensation; keep the residual
		// across encounter restarts, and consume it only on accepted steps.
		double timeCompensation = 0.0;
		int rejectionsInARow = 0;

		for (;;) {
			int completedSteps = 0;
			double lastAcceptedStep = 0.0;
			double tm = segmentStartOffset;
			if (options.cancel.stop_requested())
				return cancelledAt(tm);
			csx.setZero();
			csv.setZero();
			countedEvaluate(x, v, tm, f1);
			for (int degree = 1; degree <= 7; ++degree) {
				b[degree].setZero();
				e[degree].setZero();
			}

			bool restartEncounter = false;
			for (;;) {
				if (options.cancel.stop_requested())
					return cancelledAt(tm);
				for (auto& compensation : csb)
					compensation.setZero();
				computeIas15Beta(b, g, rc);
				const double step = tp;
				const double currentMjd = tini + tm;
				if (!std::isfinite(step))
					throw std::runtime_error("IAS15 produced a non-finite step");
				if (direction * step <= 0.0)
					throw std::runtime_error("IAS15 produced a zero or reversed step");
				if (currentMjd + step == currentMjd)
					throw std::runtime_error("IAS15 step is below Julian-date resolution");

				double correctorError = std::numeric_limits<double>::max();
				double previousError = std::numeric_limits<double>::max();

				for (int iteration = 1; iteration <= kMaxCorrectorIterations; ++iteration) {
					++report.correctorIterations;
					double maxDeltaG7 = 0.0;
					double maxAcceleration = 0.0;

					for (int j = 2; j <= 8; ++j) {
						if (options.cancel.stop_requested())
							return cancelledAt(tm);
						const double s = rc.h[j];
						predictIas15State(x, v, csx, csv, b, f1, s, step, y, z);

						countedEvaluate(y, z, tm + (s * step - timeCompensation), fj);

						for (int k = 0; k < N; ++k) {
							const double gk = (fj(k) - f1(k)) / s;
							const double oldG = g[j - 1](k);

							switch (j) {
								case 2:
									g[1](k) = gk;
									break;
								case 3:
									g[2](k) = (gk - g[1](k)) / rc.intervals[1];
									break;
								case 4:
									g[3](k) = ((gk - g[1](k)) / rc.intervals[2] - g[2](k)) / rc.intervals[3];
									break;
								case 5:
									g[4](k) =
									    (((gk - g[1](k)) / rc.intervals[4] - g[2](k)) / rc.intervals[5] - g[3](k)) /
									    rc.intervals[6];
									break;
								case 6:
									g[5](k) =
									    ((((gk - g[1](k)) / rc.intervals[7] - g[2](k)) / rc.intervals[8] - g[3](k)) /
									         rc.intervals[9] -
									     g[4](k)) /
									    rc.intervals[10];
									break;
								case 7:
									g[6](k) =
									    (((((gk - g[1](k)) / rc.intervals[11] - g[2](k)) / rc.intervals[12] - g[3](k)) /
									          rc.intervals[13] -
									      g[4](k)) /
									         rc.intervals[14] -
									     g[5](k)) /
									    rc.intervals[15];
									break;
								case 8:
									g[7](k) = ((((((gk - g[1](k)) / rc.intervals[16] - g[2](k)) / rc.intervals[17] -
									              g[3](k)) /
									                 rc.intervals[18] -
									             g[4](k)) /
									                rc.intervals[19] -
									            g[5](k)) /
									               rc.intervals[20] -
									           g[6](k)) /
									          rc.intervals[21];
									break;
								default:
									break;
							}

							const double delta = g[j - 1](k) - oldG;
							addCompensated(b[j - 1](k), csb[j - 1](k), delta);
							const int firstC = (j - 3) * (j - 2) / 2 + 1;
							for (int l = 1; l <= j - 2; ++l)
								addCompensated(b[l](k), csb[l](k), rc.c[firstC + l - 1] * delta);

							if (j == 8) {
								const double magnitude = std::abs(delta);
								const double accel = std::abs(fj(k));
								if (!std::isfinite(magnitude) || !std::isfinite(accel))
									throw std::runtime_error("IAS15 corrector produced non-finite coefficients");
								if (std::isnormal(magnitude) && magnitude > maxDeltaG7)
									maxDeltaG7 = magnitude;

								if (std::isnormal(accel) && accel > maxAcceleration)
									maxAcceleration = accel;
							}
						}
					}

					previousError = correctorError;
					correctorError = maxAcceleration > 0.0
					                     ? maxDeltaG7 / maxAcceleration
					                     : (maxDeltaG7 > 0.0 ? std::numeric_limits<double>::max() : 0.0);

					if (correctorError < kCorrectorTolerance)
						break;

					if (iteration >= kStagnationCheckIteration && previousError <= correctorError)
						break;

					if (iteration == kMaxCorrectorIterations && !fixedStep)
						++report.unconvergedEncounterSteps;
				}

				double nextStep = step;
				if (!fixedStep) {
					const double error = globalRelativeError<N>(b, fj, x, v, step);
					if (std::isnormal(error)) {
						nextStep = seventhRoot(epsilon / error) * step;
					} else {
						nextStep = step / kSafetyFactor;
					}

					if (!std::isfinite(nextStep) || nextStep == 0.0)
						nextStep = step / kSafetyFactor;

					if (std::abs(nextStep) < IAS15Integrator::kMinimumAdaptiveStepDays)
						nextStep = std::copysign(IAS15Integrator::kMinimumAdaptiveStepDays, step);

					const double ratio = std::abs(nextStep / step);

					if (ratio < kSafetyFactor) {
						++report.rejectedSteps;
						if (++rejectionsInARow > 1000)
							throw std::runtime_error("IAS15 failed to find an acceptable step");

						tp = nextStep;
						if (completedSteps > 0) {
							for (int degree = 1; degree <= 7; ++degree) {
								b[degree] = acceptedB[degree];
								e[degree] = acceptedE[degree];
							}
							predictIas15Coefficients(tp / lastAcceptedStep, b, e);
						}
						continue;
					}

					if (ratio > 1.0 / kSafetyFactor)
						nextStep = step / kSafetyFactor;
				}
				rejectionsInARow = 0;

				const RadauVector<N> stepStartX = x;
				const RadauVector<N> stepStartV = v;

				while (targetIndex < targetMjds.size()) {
					const double outputOffset = targetMjds[targetIndex] - tini;
					const double outputDelta = (outputOffset - tm) + timeCompensation;
					if (direction * (outputDelta - step) > 1.0e-10)
						break;
					if (direction * outputDelta < -1.0e-10)
						throw std::runtime_error("IAS15 skipped an output epoch");

					const double s = std::clamp(outputDelta / step, 0.0, 1.0);
					RadauVector<N> outputX;
					RadauVector<N> outputV;
					predictIas15State(stepStartX, stepStartV, csx, csv, b, f1, s, step, outputX, outputV);
					writeGlobalOutput(targetIndex, targetMjds[targetIndex], outputX, outputV);
					++targetIndex;
					report.completedOutputs = targetIndex;
				}

				for (int k = 0; k < N; ++k) {
					for (int coefficient = 7; coefficient >= 1; --coefficient)
						addCompensated(x(k), csx(k),
						               b[coefficient](k) / ((coefficient + 1) * (coefficient + 2)) * step * step);
					addCompensated(x(k), csx(k), f1(k) / 2.0 * step * step);
					addCompensated(x(k), csx(k), v(k) * step);
					for (int coefficient = 7; coefficient >= 1; --coefficient)
						addCompensated(v(k), csv(k), b[coefficient](k) / (coefficient + 1) * step);
					addCompensated(v(k), csv(k), f1(k) * step);
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
				    (observeAcceptedSteps && !observeStep(currentMjd, stepStartX, stepStartV, tini + tm, x, v))) {
					report.finalState = {Vec3D(x(0), x(1), x(2)), Vec3D(v(0), v(1), v(2))};
					report.reachedMjd = tini + tm;
					report.stoppedByObserver = true;
					return report;
				}

				if (std::abs((tm - span) - timeCompensation) < 1.0e-10) {
					if (options.cancel.stop_requested())
						return cancelledAt(tm);
					if (targetIndex != targetMjds.size())
						throw std::runtime_error("IAS15 did not emit every output");
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
						}
						++report.encounterRestarts;
						segmentStartOffset = tm;
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

				tp = fixedStep ? *fixedStep : nextStep;
				const double remaining = (span - tm) + timeCompensation;
				if (direction * (remaining - tp) < 1.0e-8)
					tp = remaining;

				for (int degree = 1; degree <= 7; ++degree) {
					acceptedB[degree] = b[degree];
					acceptedE[degree] = e[degree];
				}
				lastAcceptedStep = step;
				predictIas15Coefficients(tp / step, b, e);
			}

			if (!restartEncounter)
				throw std::runtime_error("IAS15 stopped before reaching target");
		}
	}

}  // namespace

IntegrationReport IAS15Integrator::integrateState(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
                                                  const IntegratedState& initial, double epochMjd,
                                                  std::span<const double> targetMjds,
                                                  std::span<IntegratedState> outputs, const IntegrationOptions& options,
                                                  const AcceptedStateStepObserver& observer) const
{
	RadauVector<3> x, v;
	packState(initial.R, initial.V, x, v);

	return integrateIas15<3>(
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
	    [&observer](double beginMjd, const RadauVector<3>& beginX, const RadauVector<3>& beginV, double endMjd,
	                const RadauVector<3>& endX, const RadauVector<3>& endV) {
		    return observer(AcceptedStateStep{beginMjd,
		                                      endMjd,
		                                      {positionFromRadau(beginX), velocityFromRadau(beginV)},
		                                      {positionFromRadau(endX), velocityFromRadau(endV)}});
	    },
	    static_cast<bool>(observer));
}

IntegrationReport IAS15Integrator::integrateStateWithStm(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
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

	return integrateIas15<21>(
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

IntegrationReport IAS15Integrator::integrateStateWithStmAndSens(
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

	return integrateIas15<33>(
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
