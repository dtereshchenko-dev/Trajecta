#include "fit/OrbitFitProblem.h"
#include <cmath>
#include <stdexcept>

namespace orbital_core::fit {
	void OrbitFitProblem::setIncluded(std::span<const std::uint8_t> included)
	{
		if (!included.empty() && included.size() != model_.observations().size())
			throw std::invalid_argument("Fit observation mask size mismatch");
		included_.assign(included.begin(), included.end());
	}

	bool OrbitFitProblem::isIncluded(std::size_t index) const noexcept
	{
		return included_.empty() || included_[index] != 0;
	}

	int OrbitFitProblem::activeResidualCount() const noexcept
	{
		int count = 0;
		for (std::size_t i = 0; i < model_.observations().size(); ++i)
			if (isIncluded(i))
				count += model_.observations()[i].measurementDimensions();
		return count;
	}

	void OrbitFitProblem::maskResiduals(std::vector<double>& residuals) const
	{
		for (std::size_t i = 0; i < model_.observations().size(); ++i) {
			if (!isIncluded(i)) {
				const int row = model_.residualOffset(i);
				for (int j = 0; j < model_.observations()[i].measurementDimensions(); ++j)
					residuals[row + j] = 0.0;
			}
		}
	}

	void OrbitFitProblem::maskJacobian(Eigen::MatrixXd& jacobian) const
	{
		for (std::size_t i = 0; i < model_.observations().size(); ++i)
			if (!isIncluded(i))
				jacobian.middleRows(model_.residualOffset(i), model_.observations()[i].measurementDimensions())
				    .setZero();
	}

	void OrbitFitProblem::beginPass(int parameterCount)
	{
		numParams_ = parameterCount;
		finalLinearization_ = {};
		cacheParams_.clear();
		cacheRes_.clear();
		cacheJ_.resize(0, 0);
		cacheTransforms_.clear();
		cacheUsesStm_ = false;
	}

	int OrbitFitProblem::numParams() const noexcept { return numParams_; }

	int OrbitFitProblem::numResiduals() const noexcept { return model_.numResiduals(); }

	void OrbitFitProblem::residuals(const std::vector<double>& params, std::vector<double>& outRes) const
	{
		ensureCache(params);
		outRes = cacheRes_;
		weightResiduals(outRes, cacheTransforms_);
		maskResiduals(outRes);
	}

	void OrbitFitProblem::jacobian(const std::vector<double>& params, Eigen::MatrixXd& outJ) const
	{
		ensureCache(params);
		if (cacheUsesStm_) {
			outJ = cacheJ_;
		}

		else {
			outJ = model_.computeJacobian(params);
		}
		weightJacobian(outJ, cacheTransforms_);
		maskJacobian(outJ);
	}

	bool OrbitFitProblem::acceptParams(const std::vector<double>& candidate) const
	{
		for (double x : candidate) {
			if (std::abs(x) > 1.0e15)
				return false;
		}

		constexpr double kMinSemiMajorAU = 1.0e-3;
		constexpr double kMaxEccentricity = 50.0;
		try {
			const auto orbit = model_.orbitFromStateVector(candidate);
			if (std::abs(orbit->getSemiMajorAxis()) < kMinSemiMajorAU)
				return false;

			if (orbit->getExcentricity() > kMaxEccentricity)
				return false;
		}

		catch (...) {
			return false;
		}

		return true;
	}

	double OrbitFitProblem::displayCost(const std::vector<double>& params) const
	{
		const auto residuals =
		    cacheParams_ == params && !cacheRes_.empty() ? cacheRes_ : model_.computeResiduals(params);
		if (included_.empty())
			return model_.opticalRms(residuals);
		double sum = 0.0;
		int count = 0;
		for (std::size_t i = 0; i < model_.observations().size(); ++i) {
			if (!isIncluded(i) || model_.observations()[i].radar)
				continue;
			const int row = model_.residualOffset(i);
			const double ra = residuals[row] / Constants::ARCSEC_RAD;
			const double dec = residuals[row + 1] / Constants::ARCSEC_RAD;
			sum += ra * ra + dec * dec;
			count += 2;
		}
		return count ? std::sqrt(sum / count) : 0.0;
	}

	void OrbitFitProblem::saveLinearizationFor(const std::vector<double>& params)
	{
		ensureCache(params);
		Eigen::MatrixXd jacobian = cacheUsesStm_ ? cacheJ_ : model_.computeJacobian(params);
		if (jacobian.rows() == 0 || jacobian.cols() != static_cast<Eigen::Index>(params.size()) ||
		    cacheRes_.size() != static_cast<std::size_t>(jacobian.rows())) {
			finalLinearization_ = {};
			return;
		}

		finalLinearization_.params = params;
		finalLinearization_.rawResiduals = cacheRes_;
		finalLinearization_.rawJacobian = std::move(jacobian);
		finalLinearization_.transforms = cacheTransforms_;
		finalLinearization_.valid = true;
	}

	void OrbitFitProblem::ensureCache(const std::vector<double>& params) const
	{
		if (cacheParams_ == params && (cacheUsesStm_ || !model_.usesNBody()))
			return;

		if (model_.usesNBody()) {
			auto sj = model_.computeResidualsAndStateJac(params);
			if (sj.stateJ.size() != 0) {
				cacheRes_ = std::move(sj.residuals);
				cacheJ_ = std::move(sj.stateJ);
				cacheTransforms_ = std::move(sj.transforms);
				cacheUsesStm_ = true;
			}

			else {
				cacheRes_ = std::move(sj.residuals);
				cacheJ_.resize(0, 0);
				cacheTransforms_.clear();
				cacheUsesStm_ = false;
			}
		}

		else {
			cacheRes_ = model_.computeResiduals(params);
			cacheJ_.resize(0, 0);
			cacheTransforms_.clear();
			cacheUsesStm_ = false;
		}
		cacheParams_ = params;
	}

	obs::AstrometricResidualPair OrbitFitProblem::normalizedResidualPair(
	    int i, const std::vector<double>& res, std::span<const obs::AstrometricResidualTransform> transforms) const
	{
		if (i < 0 || static_cast<std::size_t>(i) >= model_.observations().size())
			return {};

		const int row = model_.residualOffset(i);
		if (res.size() < static_cast<std::size_t>(row + 2))
			return {};

		if (transforms.size() == model_.observations().size() && transforms[static_cast<std::size_t>(i)].valid) {
			return orbital_core::obs::applyAstrometricResidualTransform(res[row], res[row + 1],
			                                                            transforms[static_cast<std::size_t>(i)]);
		}

		return {};
	}

	void OrbitFitProblem::weightResiduals(std::vector<double>& res,
	                                      std::span<const obs::AstrometricResidualTransform> transforms) const
	{
		const int N = static_cast<int>(model_.observations().size());
		if (!model_.usesWeights())
			return;

		for (int i = 0; i < N; ++i) {
			int row = model_.residualOffset(i);
			const auto& observation = model_.observations()[i];
			if (observation.radar) {
				for (const auto& measurement : {observation.radar->delaySeconds, observation.radar->dopplerHz})
					if (measurement)
						res[row++] /= measurement->sigma;
				continue;
			}
			const auto pair = normalizedResidualPair(i, res, transforms);
			if (pair.valid) {
				res[row] = pair.first;
				res[row + 1] = pair.second;
				continue;
			}

			const double wr = 1.0 / model_.observations()[i].sigmaRa;
			const double wd = 1.0 / model_.observations()[i].sigmaDec;
			res[row] *= wr;
			res[row + 1] *= wd;
		}
	}

	double OrbitFitProblem::sumSquares(const std::vector<double>& res) const
	{
		double sumSq = 0.0;
		for (std::size_t i = 0; i < model_.observations().size(); ++i) {
			const int row = model_.residualOffset(i);
			if (model_.observations()[i].measurementDimensions() == 2)
				sumSq += res[row] * res[row] + res[row + 1] * res[row + 1];
			else
				sumSq += res[row] * res[row];
		}
		return sumSq;
	}

	void OrbitFitProblem::weightJacobian(Eigen::MatrixXd& J,
	                                     std::span<const obs::AstrometricResidualTransform> transforms) const
	{
		const int N = static_cast<int>(model_.observations().size());
		if (!model_.usesWeights())
			return;

		for (int i = 0; i < N; ++i) {
			int row = model_.residualOffset(i);
			const auto& observation = model_.observations()[i];
			if (observation.radar) {
				for (const auto& measurement : {observation.radar->delaySeconds, observation.radar->dopplerHz})
					if (measurement)
						J.row(row++) /= measurement->sigma;
				continue;
			}
			if (transforms.size() == model_.observations().size() && transforms[static_cast<std::size_t>(i)].valid) {
				const auto& tr = transforms[static_cast<std::size_t>(i)];
				const Eigen::RowVectorXd row0 = J.row(row);
				const Eigen::RowVectorXd row1 = J.row(row + 1);
				J.row(row) = tr.m00 * row0 + tr.m01 * row1;
				J.row(row + 1) = tr.m10 * row0 + tr.m11 * row1;
				continue;
			}

			const double wr = 1.0 / model_.observations()[i].sigmaRa;
			const double wd = 1.0 / model_.observations()[i].sigmaDec;
			J.row(row) *= wr;
			J.row(row + 1) *= wd;
		}
	}

	WeightedObservationLinearization OrbitFitProblem::weightedLinearization(const std::vector<double>& params) const
	{
		WeightedObservationLinearization out;

		std::vector<double> residuals;
		Eigen::MatrixXd jacobian;
		std::vector<orbital_core::obs::AstrometricResidualTransform> transforms;

		if (finalLinearization_.matches(params)) {
			residuals = finalLinearization_.rawResiduals;
			jacobian = finalLinearization_.rawJacobian;
			transforms = finalLinearization_.transforms;
		} else if (model_.usesNBody()) {
			auto sj = model_.computeResidualsAndStateJac(params);
			residuals = std::move(sj.residuals);
			jacobian = std::move(sj.stateJ);
			transforms = std::move(sj.transforms);
		}

		if (residuals.empty())
			residuals = model_.computeResiduals(params);

		if (jacobian.size() == 0) {
			jacobian = model_.computeJacobian(params);
		}

		if (residuals.size() != model_.numResiduals() || jacobian.rows() != model_.numResiduals() ||
		    jacobian.cols() != static_cast<Eigen::Index>(params.size())) {
			return out;
		}

		weightResiduals(residuals, transforms);
		const double sumSq = sumSquares(residuals);

		weightJacobian(jacobian, transforms);

		jacobian.col(3) *= 1.0 / Constants::GAUSS_K;
		jacobian.col(4) *= 1.0 / Constants::GAUSS_K;
		jacobian.col(5) *= 1.0 / Constants::GAUSS_K;

		out.residuals = std::move(residuals);
		out.jacobian = std::move(jacobian);
		out.sumSquares = sumSq;
		out.valid = true;
		return out;
	}

	double OrbitFitProblem::weightedChiSquare(const std::vector<double>& params) const
	{
		if (model_.usesNBody()) {
			auto sj = model_.computeResidualsAndStateJac(params);
			if (!sj.residuals.empty()) {
				weightResiduals(sj.residuals, sj.transforms);
				return sumSquares(sj.residuals);
			}
		}

		return model_.computeCostWeighted(params);
	}

}  // namespace orbital_core::fit
