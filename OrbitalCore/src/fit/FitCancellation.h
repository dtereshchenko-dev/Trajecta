#pragma once

#include "orbital_core/FitRunControl.h"

#include <exception>
#include <stop_token>

namespace orbital_core::fit {

	// Internal unwinding only; all fit entry points return empty cancelled results.
	// Cancellation is cooperative between evaluations, not inside the integrator.
	class FitCancelled final : public std::exception {
	   public:
		const char* what() const noexcept override { return "Fit cancelled"; }
	};

	inline void checkCancellation(std::stop_token cancel)
	{
		if (cancel.stop_requested())
			throw FitCancelled{};
	}

	// A callback may request stop; do not enter the announced stage afterwards.
	inline void reportProgress(const FitRunControl& control, const std::string& message)
	{
		checkCancellation(control.cancel);
		if (control.progress)
			control.progress(message);
		checkCancellation(control.cancel);
	}

}  // namespace orbital_core::fit
