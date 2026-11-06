#pragma once

#include <cmath>
#include <optional>
#include <string>

namespace orbital_core::obs {

	struct RadarMeasurement {
		double value{};
		double sigma{};
		bool operator==(const RadarMeasurement&) const = default;
	};

	enum class RadarReturn { CenterOfMass, Surface };

	// Reception epoch is the observation UTC. Delay is round-trip seconds;
	// Doppler is received minus transmitted frequency, Hz (approaching is positive).
	// Missing channels are absent, never zero-valued synthetic measurements.
	struct RadarObservation {
		std::optional<RadarMeasurement> delaySeconds;
		std::optional<RadarMeasurement> dopplerHz;
		double frequencyHz{};
		// RWO gives range rate, not transmitter frequency. For such imports,
		// Doppler and its sigma are represented at a reference frequency of 1 Hz.
		// This is a unit conversion, not an inferred transmitter frequency.
		bool frequencyIsReference{false};
		std::string transmitter;
		RadarReturn returnType{RadarReturn::CenterOfMass};
		// Surface returns only; preparation assumes 1 km when absent.
		std::optional<double> surfaceRadiusKm;

		[[nodiscard]] int dimensions() const noexcept
		{
			return static_cast<int>(delaySeconds.has_value()) + static_cast<int>(dopplerHz.has_value());
		}

		[[nodiscard]] bool valid() const noexcept
		{
			const auto validChannel = [](const auto& channel) {
				return !channel ||
				       (std::isfinite(channel->value) && std::isfinite(channel->sigma) && channel->sigma > 0.0);
			};
			return dimensions() > 0 && !transmitter.empty() && validChannel(delaySeconds) && validChannel(dopplerHz) &&
			       (!delaySeconds || delaySeconds->value > 0.0) && std::isfinite(frequencyHz) && frequencyHz > 0.0 &&
			       (!surfaceRadiusKm || (std::isfinite(*surfaceRadiusKm) && *surfaceRadiusKm > 0.0));
		}

		bool operator==(const RadarObservation&) const = default;
	};

}  // namespace orbital_core::obs
