#include "hazard/lov/LovTraceDetail.h"
#include <optional>

#include "hazard/HazardAudit.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>

namespace orbital_core::hazard::lov_detail {
	namespace {
		LovMinimumResult finishMinimum(const char* method, const TraceRecord& start, std::optional<TraceRecord> record,
		                               bool converged, bool noTargetPlane, int iterations, bool reclassified = false)
		{
			if (record)
				record->riskCandidate = !noTargetPlane;
			LOG_INFO << std::setprecision(17) << "[LineOfVariation] minimum-result"
			         << " method=" << method << " sample=" << start.sampleIndex << " startSigma=" << start.trace.sigma
			         << " iterations=" << iterations << " noTp=" << (noTargetPlane ? "yes" : "no")
			         << " converged=" << (converged ? "yes" : "no") << " reclassified=" << (reclassified ? "yes" : "no")
			         << " hasRecord=" << (record ? "yes" : "no") << " sigma=" << (record ? record->trace.sigma : 0.0)
			         << " jd=" << (record ? record->event.jd : 0.0)
			         << " b=" << (record ? record->trace.distanceBodyRadii : 0.0)
			         << " stretch=" << (record ? record->trace.stretchBodyRadii : 0.0)
			         << " p=" << (record ? record->probability.probability : 0.0);
			return {std::move(record), converged, noTargetPlane, reclassified};
		}

		void logBracketInput(const TraceRecord& left, const TraceRecord& right)
		{
			LOG_INFO << std::setprecision(17) << "[LineOfVariation] minimum-bracket"
			         << " leftSample=" << left.sampleIndex << " rightSample=" << right.sampleIndex
			         << " leftSigma=" << left.trace.sigma << " rightSigma=" << right.trace.sigma
			         << " leftJd=" << left.event.jd << " rightJd=" << right.event.jd
			         << " leftB=" << left.trace.distanceBodyRadii << " rightB=" << right.trace.distanceBodyRadii
			         << " leftDerivative=" << left.trace.dDistance2DSigmaBodyRadii2
			         << " rightDerivative=" << right.trace.dDistance2DSigmaBodyRadii2;
		}
	}  // namespace

	LovMinimumResult refineBracketedMinimum(const TraceRecord& left, const TraceRecord& right,
	                                        const LovEncounterEvaluator& evaluate)
	{
		constexpr int kMaxIterations = 30;
		constexpr double kDistanceTolerance = 2.3e-3;
		constexpr double kProbabilityTolerance = 1.0e-10;

		double a = left.trace.sigma;
		double b = right.trace.sigma;
		double fa = left.trace.dDistance2DSigmaBodyRadii2;
		double fb = right.trace.dDistance2DSigmaBodyRadii2;
		double da = left.trace.distanceBodyRadii;
		double db = right.trace.distanceBodyRadii;
		double ta = left.event.jd;
		double tb = right.event.jd;
		TraceRecord current = left;
		int iterations = 0;
		double previousDistance = da;

		if (!(a < b) || !(fa < 0.0) || !(fb > 0.0)) {
			return finishMinimum("falsi", left, {}, false, true, 0);
		}

		double w = (fb * a - fa * b) / (fb - fa);
		double previousTrialF = fa;
		double previousPreviousTrialF = std::numeric_limits<double>::quiet_NaN();

		logBracketInput(left, right);
		for (int iter = 0; iter < kMaxIterations; ++iter) {
			iterations = iter + 1;
			if (!(w > a && w < b))
				w = 0.5 * (a + b);

			std::optional<TraceRecord> mid = evaluate(left.sampleIndex, w, ta, tb);
			if (!mid) {
				return finishMinimum("falsi", left, current, false, true, iterations);
			}

			const double fw = mid->trace.dDistance2DSigmaBodyRadii2;
			current = *mid;
			double effectiveFa = fa;
			double effectiveFb = fb;
			if (fa * fw <= 0.0) {
				b = w;
				fb = fw;
				db = mid->trace.distanceBodyRadii;
				tb = mid->event.jd;

				if (iter > 0 && previousTrialF * fw > 0.0) {
					effectiveFa = fa * (previousPreviousTrialF * previousTrialF > 0.0 ? 0.25 : 0.5);
				}
			} else {
				a = w;
				fa = fw;
				da = mid->trace.distanceBodyRadii;
				ta = mid->event.jd;

				if (iter > 0 && previousTrialF * fw > 0.0) {
					effectiveFb = fb * (previousPreviousTrialF * previousTrialF > 0.0 ? 0.25 : 0.5);
				}
			}

			const double deltaSigma = std::abs(b - a);
			const double distanceControlByDerivative = std::abs(fw) * deltaSigma / std::max(da + db, 1.0e-30);
			const double distanceControlByStep = std::abs(mid->trace.distanceBodyRadii - previousDistance);
			const double distanceControl =
			    iter == 0 ? distanceControlByDerivative : std::min(distanceControlByDerivative, distanceControlByStep);
			const double probabilityControl = deltaSigma / 6.0;
			previousDistance = mid->trace.distanceBodyRadii;

			if (distanceControl < kDistanceTolerance || probabilityControl < kProbabilityTolerance ||
			    current.trace.distanceBodyRadii < current.trace.impactCrossSectionBodyRadii) {
				return finishMinimum("falsi", left, current, true, false, iterations);
			}

			previousPreviousTrialF = previousTrialF;
			previousTrialF = fw;
			const double denom = effectiveFb - effectiveFa;
			if (std::abs(denom) < 1.0e-30)
				break;
			w = (effectiveFb * a - effectiveFa * b) / denom;
		}

		return finishMinimum("falsi", left, current, current.trace.distanceBodyRadii < 1.0, false, iterations);
	}

	LovMinimumResult refineEndpointMinimum(const TraceRecord& start, double sigmaStepLimit,
	                                       const LovEncounterEvaluator& evaluate)
	{
		constexpr int kMaxIterations = 8;
		constexpr double kDistanceTolerance = 2.3e-3;
		constexpr double kProbabilityTolerance = 1.0e-10;

		if (!(sigmaStepLimit >= 0.0) || !std::isfinite(sigmaStepLimit) || !start.trace.valid) {
			return finishMinimum("achilles", start, {}, false, true, 0);
		}

		TraceRecord current = start;
		double stepLimit = sigmaStepLimit;
		bool noTargetPlane = false;
		bool reclassified = false;
		int iterations = 0;
		audit::Scope auditScope(start.sampleIndex, start.trace.sigma);
		audit::write("endpoint-begin", [&](auto& log) {
			log << " jd=" << start.event.jd << " b=" << start.trace.distanceBodyRadii
			    << " derivative=" << start.trace.dDistance2DSigmaBodyRadii2
			    << " stretchLov=" << start.trace.stretchLovBodyRadiiPerSigma << " stepLimit=" << stepLimit;
		});

		for (int iter = 0; iter < kMaxIterations; ++iter) {
			iterations = iter + 1;
			const double stretch = current.trace.stretchLovBodyRadiiPerSigma;
			if (!(stretch > 0.0) || !std::isfinite(stretch)) {
				return finishMinimum("achilles", start, {}, false, true, iterations, reclassified);
			}

			const double curvature = 2.0 * stretch * stretch;
			double step = -current.trace.dDistance2DSigmaBodyRadii2 / curvature;
			step = std::clamp(step, -stepLimit, stepLimit);
			audit::write("endpoint-step", [&](auto& log) {
				log << " iteration=" << iter + 1 << " currentSigma=" << current.trace.sigma
				    << " trialSigma=" << current.trace.sigma + step << " step=" << step << " stepLimit=" << stepLimit
				    << " referenceJd=" << current.event.jd;
			});

			std::optional<TraceRecord> next =
			    evaluate(current.sampleIndex, current.trace.sigma + step, current.event.jd, current.event.jd);
			noTargetPlane = !next.has_value();
			if (noTargetPlane) {
				audit::write("endpoint-result",
				             [&](auto& log) { log << " iteration=" << iter + 1 << " status=no-tp"; });
				stepLimit = std::min(0.5 * stepLimit, 0.5 * std::abs(step));
				continue;
			}

			next->riskCandidate = true;
			audit::write("endpoint-result", [&](auto& log) {
				log << " iteration=" << iter + 1 << " status=ok"
				    << " trialSigma=" << next->trace.sigma << " jd=" << next->event.jd
				    << " b=" << next->trace.distanceBodyRadii
				    << " impactRadius=" << next->trace.impactCrossSectionBodyRadii
				    << " derivative=" << next->trace.dDistance2DSigmaBodyRadii2
				    << " stretchLov=" << next->trace.stretchLovBodyRadiiPerSigma;
			});

			const double d0 = current.trace.distanceBodyRadii;
			const double d1 = next->trace.distanceBodyRadii;
			const double distanceControlByDerivative =
			    std::abs(next->trace.dDistance2DSigmaBodyRadii2) * std::abs(step) / std::max(d0 + d1, 1.0e-30);
			const double distanceControlByStep = std::abs(d1 - d0);
			const double distanceControl =
			    iter == 0 ? distanceControlByDerivative : std::min(distanceControlByDerivative, distanceControlByStep);
			const double probabilityControl = std::abs(step) / 6.0;

			if (distanceControl < kDistanceTolerance || distanceControl < kDistanceTolerance * d1 ||
			    probabilityControl < kProbabilityTolerance || d1 < next->trace.impactCrossSectionBodyRadii) {
				return finishMinimum("achilles", start, std::move(next), true, false, iterations, reclassified);
			}

			const TraceRecord* lower = &current;
			const TraceRecord* upper = &*next;
			if (lower->trace.sigma > upper->trace.sigma)
				std::swap(lower, upper);

			if (lower->trace.dDistance2DSigmaBodyRadii2 < 0.0 && upper->trace.dDistance2DSigmaBodyRadii2 > 0.0) {
				reclassified = true;
				if (classifyLovPair(*lower, *upper) == LovPairType::InterruptedMinimum) {
					LovMinimumResult bracketed = refineBracketedMinimum(*lower, *upper, evaluate);
					return finishMinimum("achilles-falsi", start, std::move(bracketed.record), bracketed.converged,
					                     bracketed.noTargetPlane, iterations, reclassified);
				}
			}

			current = std::move(*next);
			stepLimit *= 0.5;
		}

		return finishMinimum("achilles", start, std::move(current), false, noTargetPlane, iterations, reclassified);
	}

}  // namespace orbital_core::hazard::lov_detail
