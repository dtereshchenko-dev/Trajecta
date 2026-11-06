#include "hazard/lov/LovTraceDetail.h"
#include <cmath>

#include "orbital_core/Constants.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/nongrav/NonGravProfile.h"

#include <limits>
#include <vector>

namespace orbital_core::hazard::lov_detail {

	static Eigen::Matrix<double, 6, 1> nominalState(const Orbit& orbit)
	{
		const Vec3D r = orbit.getR();
		const Vec3D v = orbit.getV() * Constants::GAUSS_K;
		Eigen::Matrix<double, 6, 1> x;
		x << r.x, r.y, r.z, v.x, v.y, v.z;
		return x;
	}

	VecX nominalParameters(const Orbit& orbit, const NonGravModel& ng)
	{
		const int nExtra = ng.numFitParams();
		VecX x(6 + nExtra);
		x.head<6>() = nominalState(orbit);

		if (nExtra == 0)
			return x;

		for (int k = 0; k < nExtra; ++k)
			x(6 + k) = ng.fitParameterValue(k);
		return x;
	}

	NonGravModel nonGravFromParameters(const NonGravModel& nominal, const VecX& params)
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

	LovImpactProbabilityResult probabilityForTrace(const TargetPlaneTracePoint& trace)
	{
		LovImpactProbabilityInput input;
		input.tracePoint = trace;
		input.impactCrossSectionBodyRadii = trace.impactCrossSectionBodyRadii;
		return lovImpactProbability(input);
	}

	bool validConfig(const LineOfVariationImpactConfig& cfg) noexcept
	{
		return cfg.closeApproach.horizonYears > 0.0 && cfg.closeApproach.lookbackYears >= 0.0 &&
		       cfg.closeApproach.thresholdAu > 0.0 && cfg.closeApproach.encounterRadiusAu > 0.0 &&
		       cfg.closeApproach.sampleStepDay > 0.0 && cfg.sigmaMax > 0.0 && cfg.sigmaStepMax > 0.0 &&
		       (!cfg.probabilitySampling || cfg.probabilityStep > 0.0) && cfg.encounterStepPoints >= 0 &&
		       cfg.maxPositiveSamples > 0;
	}

	double nominalObservationChiSquareRaw(const TraceEvaluationContext& ctx)
	{
		const ImpactRiskRmsScaling& scaling = ctx.cfg.impactRefinement.rmsScaling;
		if (!scaling.linearization)
			return std::numeric_limits<double>::quiet_NaN();

		const VecX xNom = nominalParameters(ctx.nominalOrbit, ctx.nominalNonGrav);
		const auto lin = scaling.linearization(xNom);
		if (!lin.valid)
			return std::numeric_limits<double>::quiet_NaN();

		return lin.sumSquares;
	}

	double nominalObservationChiSquare(const TraceEvaluationContext& ctx)
	{
		const double raw = nominalObservationChiSquareRaw(ctx);
		if (!std::isfinite(raw))
			return raw;

		const ImpactRiskRmsScaling& scaling = ctx.cfg.impactRefinement.rmsScaling;

		const double scale =
		    std::isfinite(scaling.objectiveScale) && scaling.objectiveScale > 0.0 ? scaling.objectiveScale : 1.0;
		return raw / scale;
	}

}  // namespace orbital_core::hazard::lov_detail
