#include "PerturberSelector.h"
#include <span>
#include <cstdint>

#include <algorithm>
#include <cmath>

namespace Calculations {
	namespace {

		constexpr double kFoEarthPerturberRangeAu = 0.2;

		constexpr auto kMercury = orbital_core::planetaryBodyMask(EphemerisPlanet::Mercury);
		constexpr auto kVenus = orbital_core::planetaryBodyMask(EphemerisPlanet::Venus);
		constexpr auto kEarth = orbital_core::planetaryBodyMask(EphemerisPlanet::Earth);
		constexpr auto kMars = orbital_core::planetaryBodyMask(EphemerisPlanet::Mars);
		constexpr auto kJupiter = orbital_core::planetaryBodyMask(EphemerisPlanet::Jupiter);
		constexpr auto kSaturn = orbital_core::planetaryBodyMask(EphemerisPlanet::Saturn);
		constexpr auto kUranus = orbital_core::planetaryBodyMask(EphemerisPlanet::Uranus);
		constexpr auto kNeptune = orbital_core::planetaryBodyMask(EphemerisPlanet::Neptune);
		constexpr auto kPluto = orbital_core::planetaryBodyMask(EphemerisPlanet::Pluto);

		[[nodiscard]] uint16_t scheduledMaskForArc(double arcDays) noexcept
		{
			uint16_t mask = 0;
			if (arcDays > 10.0)
				mask |= kJupiter;
			if (arcDays > 40.0)
				mask |= kSaturn;
			if (arcDays > 200.0)
				mask |= kVenus | kEarth | kMars | kUranus;
			if (arcDays > 365.25 * 5.0)
				mask |= kMercury | kNeptune | kPluto;
			return mask;
		}

	}  // namespace

	PerturberSelection ArcScheduledPerturberSelector::select(double arcDays,
	                                                         uint16_t carriedAutomaticMask) const noexcept
	{
		PerturberSelection out;

		if (!m_config.nbodyEnabled)
			return out;

		if (m_config.mode == PerturberSelectionMode::ManualExact) {
			out.effectiveMask = m_config.userMask;
			return out;
		}

		out.scheduledMask = scheduledMaskForArc(arcDays);
		out.automaticMask = carriedAutomaticMask;
		out.alwaysIncludedMask = m_config.alwaysIncludedMask;

		const uint16_t requested =
		    static_cast<uint16_t>(out.scheduledMask | out.automaticMask | out.alwaysIncludedMask);
		out.effectiveMask = static_cast<uint16_t>(requested & m_config.userMask);
		out.blockedMask = static_cast<uint16_t>(requested & ~m_config.userMask);
		return out;
	}

	uint16_t ArcScheduledPerturberSelector::detectAutomaticMask(
	    std::span<const Vec3D> positions,
	    std::span<const orbital_core::obs::PreparedObservation> prepared) const noexcept
	{
		if (!m_config.nbodyEnabled || m_config.mode != PerturberSelectionMode::AutoByArcLength)
			return 0;

		const std::size_t count = std::min(positions.size(), prepared.size());
		for (std::size_t i = 0; i < count; ++i) {
			const double geocentricDistance = (positions[i] - prepared[i].earthR).magnitude();
			if (std::isfinite(geocentricDistance) && geocentricDistance < kFoEarthPerturberRangeAu) {
				return static_cast<uint16_t>(kEarth | MOON_BODY_MASK);
			}
		}
		return 0;
	}

}  // namespace Calculations
