#pragma once

#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include <Eigen/Core>

namespace orbital_core::hazard::mc {

	// Freezes measurement-noise factors using the same full-pool preparation,
	// final NG model and optical timing/correlation transform as the accepted fit.
	// The nominal covariance scale is not a measurement-noise multiplier.
	class ObservationSampler {
	   public:
		ObservationSampler(const FitPipelineResult& source, std::stop_token cancel);
		[[nodiscard]] FitResult refit(int index, std::uint64_t seed, std::stop_token cancel) const;

	   private:
		const FitPipelineResult& source_;
		std::vector<std::uint8_t> included_;
		std::vector<Eigen::Matrix<double, 2, 2, Eigen::DontAlign>> factors_;
		std::vector<double> cosDec_;
	};

}  // namespace orbital_core::hazard::mc
