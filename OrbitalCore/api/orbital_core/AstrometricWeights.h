#pragma once

#include "orbital_core/Constants.h"

namespace orbital_core {

	enum class AstrometricWeightProfile : int {
		Vfcc17 = 1,
		Uniform = 2,
	};

	inline constexpr double kDefaultTimingSigmaSeconds = 0.0;

	struct ObservationWeightingOptions {
		// Prefer the observer's own ADES rmsRA/rmsDec over the table value
		// where the batch supplies one.
		bool useObserverSigmas{false};

		// Timing-uncertainty term, in seconds.  A negative value disables it
		// entirely (the historical behaviour).  Zero means "use the ADES
		// rmsTime where the observer reports one and nothing elsewhere"; a
		// positive value is additionally applied to every observation that
		// reports none.
		//
		// The term itself is the standard one: a timing error displaces the
		// object along its apparent motion, so sigma_t enters as
		//   cov += sigma_t^2 * (v v^T),   v = (d(RA cos DEC)/dt, d(DEC)/dt)
		// which our error-ellipse form already expresses as a tilt.
		double timingSigmaDefaultSeconds{-1.0};

		[[nodiscard]] constexpr bool timingEnabled() const noexcept { return timingSigmaDefaultSeconds >= 0.0; }
	};

	struct ObservationUncertainty {
		double sigmaRaRad{Constants::ARCSEC_RAD};
		double sigmaDecRad{Constants::ARCSEC_RAD};
		bool forceWeightRa{false};
		bool forceWeightDec{false};
		bool resolved{false};
	};

}  // namespace orbital_core
