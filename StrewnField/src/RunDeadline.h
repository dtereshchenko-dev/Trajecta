#pragma once

#include "strewn_core/Types.h"
#include <cmath>
#include <stdexcept>

namespace strewn::detail {

	using Deadline = std::chrono::steady_clock::time_point;

	inline void validateRunTimeout(double seconds)
	{
		if (!std::isfinite(seconds) || seconds <= 0 || seconds > 86400)
			throw std::invalid_argument("run_timeout_s must be positive and <= 86400 seconds");
	}

	inline Deadline runDeadline(double seconds)
	{
		validateRunTimeout(seconds);
		return std::chrono::steady_clock::now() +
		       std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(seconds));
	}

	inline std::optional<FlightEvent> interruption(std::stop_token cancel, Deadline deadline)
	{
		if (cancel.stop_requested())
			return FlightEvent::Cancelled;
		if (deadline != Deadline::max() && std::chrono::steady_clock::now() >= deadline)
			return FlightEvent::WallTimeout;
		return std::nullopt;
	}

	class RunInterrupted final : public std::runtime_error {
	   public:
		explicit RunInterrupted(FlightEvent event)
		    : std::runtime_error(event == FlightEvent::Cancelled ? "Calculation cancelled"
		                                                         : "Realization wall-clock timeout"),
		      m_status(event == FlightEvent::Cancelled ? RunStatus::Cancelled : RunStatus::WallTimeout)
		{
		}

		RunStatus status() const noexcept { return m_status; }

	   private:
		RunStatus m_status;
	};

	inline void checkInterruption(std::stop_token cancel, Deadline deadline)
	{
		if (const auto event = interruption(cancel, deadline))
			throw RunInterrupted(*event);
	}

}  // namespace strewn::detail
