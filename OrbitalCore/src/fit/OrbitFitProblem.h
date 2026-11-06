#pragma once

#include "fit/OrbitFitLinearization.h"
#include "fit/OrbitFitModel.h"
#include "orbital_core/optimizers/ILSQProblem.h"
#include <span>
#include <cstdint>

namespace orbital_core::fit {

	// Solver adapter and its caches. beginPass invalidates all pass-local state;
	// the accepted snapshot is shared with covariance and residual diagnostics.
	class OrbitFitProblem final : public optimizers::ILSQProblem {
	   public:
		explicit OrbitFitProblem(const OrbitFitModel& model) : model_(model) {}

		void beginPass(int parameterCount);
		// Keep raw rows for diagnostics; the solver sees zeroed rows for excluded observations.
		void setIncluded(std::span<const std::uint8_t> included);
		[[nodiscard]] bool isIncluded(std::size_t index) const noexcept;

		[[nodiscard]] std::span<const std::uint8_t> included() const noexcept { return included_; }

		[[nodiscard]] int activeResidualCount() const noexcept;
		void maskResiduals(std::vector<double>& residuals) const;
		void maskJacobian(Eigen::MatrixXd& jacobian) const;
		int numParams() const noexcept override;
		int numResiduals() const noexcept override;
		void residuals(const std::vector<double>& params, std::vector<double>& out) const override;
		void jacobian(const std::vector<double>& params, Eigen::MatrixXd& out) const override;
		bool acceptParams(const std::vector<double>& candidate) const override;
		double displayCost(const std::vector<double>& params) const override;
		void saveLinearizationFor(const std::vector<double>& params);

		[[nodiscard]] const LinearizationSnapshot& snapshot() const noexcept { return finalLinearization_; }

		[[nodiscard]] WeightedObservationLinearization weightedLinearization(const std::vector<double>& params) const;
		[[nodiscard]] double weightedChiSquare(const std::vector<double>& params) const;
		void weightJacobian(Eigen::MatrixXd& jacobian,
		                    std::span<const obs::AstrometricResidualTransform> transforms) const;
		void weightResiduals(std::vector<double>& residuals,
		                     std::span<const obs::AstrometricResidualTransform> transforms) const;
		double sumSquares(const std::vector<double>& weightedResiduals) const;

	   private:
		void ensureCache(const std::vector<double>& params) const;
		obs::AstrometricResidualPair normalizedResidualPair(
		    int index, const std::vector<double>& residuals,
		    std::span<const obs::AstrometricResidualTransform> transforms) const;

		const OrbitFitModel& model_;
		int numParams_{6};
		std::vector<std::uint8_t> included_;
		LinearizationSnapshot finalLinearization_;
		mutable std::vector<double> cacheParams_;
		mutable std::vector<double> cacheRes_;
		mutable std::vector<obs::AstrometricResidualTransform> cacheTransforms_;
		mutable Eigen::MatrixXd cacheJ_;
		mutable bool cacheUsesStm_{false};
	};

}  // namespace orbital_core::fit
