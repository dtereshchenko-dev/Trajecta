#pragma once

#include "NBodyPropagator.h"
#include "obs/PreparedObservation.h"
#include "orbital_core/Vec3D.h"

#include <cstdint>
#include <span>

namespace Calculations {

	enum class PerturberSelectionMode : std::uint8_t {
		ManualExact,
		AutoByArcLength,
	};

	struct PerturberSelectionConfig {
		bool nbodyEnabled{true};
		uint16_t userMask{static_cast<uint16_t>(PLANET_BODY_MASK | MOON_BODY_MASK)};
		uint16_t alwaysIncludedMask{0};
		PerturberSelectionMode mode{PerturberSelectionMode::AutoByArcLength};
	};

	struct PerturberSelection {
		uint16_t scheduledMask{0};
		uint16_t automaticMask{0};
		uint16_t alwaysIncludedMask{0};
		uint16_t blockedMask{0};
		uint16_t effectiveMask{0};
	};

	class ArcScheduledPerturberSelector {
	   public:
		explicit constexpr ArcScheduledPerturberSelector(PerturberSelectionConfig config) noexcept : m_config(config) {}

		[[nodiscard]] PerturberSelection select(double arcDays, uint16_t carriedAutomaticMask) const noexcept;

		[[nodiscard]] uint16_t detectAutomaticMask(
		    std::span<const Vec3D> positions,
		    std::span<const orbital_core::obs::PreparedObservation> prepared) const noexcept;

		[[nodiscard]] constexpr uint16_t userMask() const noexcept { return m_config.userMask; }

		[[nodiscard]] constexpr bool nbodyEnabled() const noexcept { return m_config.nbodyEnabled; }

	   private:
		PerturberSelectionConfig m_config;
	};

}  // namespace Calculations
