#include "hazard/LineOfVariation.h"
#include <vector>

#include "orbital_core/Constants.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::hazard {

	LovWeakDirection weakDirectionFromCovariance(const Eigen::MatrixXd& covariance)
	{
		LovWeakDirection out;

		if (covariance.rows() == 0 || covariance.rows() != covariance.cols() || !covariance.array().isFinite().all())
			return out;

		const Eigen::MatrixXd sym = 0.5 * (covariance + covariance.transpose());

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sym);
		if (solver.info() != Eigen::Success)
			return out;

		const Eigen::VectorXd values = solver.eigenvalues();
		const Eigen::MatrixXd vectors = solver.eigenvectors();
		const Eigen::Index idx = values.size() - 1;
		const double lambda = values(idx);
		if (!std::isfinite(lambda) || lambda <= 0.0)
			return out;

		Eigen::VectorXd axis = vectors.col(idx);
		for (Eigen::Index i = 0; i < axis.size(); ++i) {
			if (std::abs(axis(i)) > 1.0e-15) {
				if (axis(i) < 0.0)
					axis = -axis;
				break;
			}
		}

		out.parameterStepPerSigma = axis * std::sqrt(lambda);
		out.covarianceEigenvalue = lambda;
		out.valid = out.parameterStepPerSigma.array().isFinite().all();
		return out;
	}

	std::vector<double> sampleLovSigmas(const LovSamplingConfig& config)
	{
		std::vector<double> positive;
		positive.push_back(0.0);

		if (!std::isfinite(config.sigmaMax) || config.sigmaMax <= 0.0 || !std::isfinite(config.sigmaStepMax) ||
		    config.sigmaStepMax <= 0.0 || config.maxPositiveSamples <= 0)
			return positive;

		if (config.probabilitySampling) {
			if (!std::isfinite(config.probabilityStep) || config.probabilityStep <= 0.0 ||
			    !std::isfinite(config.diskRadiusBodyRadii) || config.diskRadiusBodyRadii <= 0.0)
				return positive;

			constexpr double kInvSqrtTwoPi = 0.398942280401432677939946059934381868;

			const double halfDiskProbability = 0.5 * config.diskRadiusBodyRadii * config.probabilityStep;

			double sigma = 0.0;
			double deltaSigma = halfDiskProbability / kInvSqrtTwoPi;

			for (int i = 0; i < config.maxPositiveSamples; ++i) {
				sigma += deltaSigma;
				if (sigma > config.sigmaMax)
					break;

				positive.push_back(sigma);

				const double pdf = kInvSqrtTwoPi * std::exp(-0.5 * sigma * sigma);
				if (pdf <= std::numeric_limits<double>::min())
					break;

				deltaSigma = std::min(halfDiskProbability / pdf, config.sigmaStepMax);
			}
		} else {
			for (int i = 0; i < config.maxPositiveSamples; ++i) {
				const double sigma = (static_cast<double>(i) + 1.0) * config.sigmaStepMax;
				if (sigma > config.sigmaMax)
					break;
				positive.push_back(sigma);
			}
		}

		std::vector<double> out;
		out.reserve(positive.size() * 2 - 1);
		for (auto it = positive.rbegin(); it != positive.rend(); ++it) {
			if (*it != 0.0)
				out.push_back(-*it);
		}
		out.insert(out.end(), positive.begin(), positive.end());
		return out;
	}

}  // namespace orbital_core::hazard
