#include "EarthOrientation.h"
#include <filesystem>
#include <optional>
#include "EopData.h"
#include "orbital_core/Types.h"
#include "orbital_core/RuntimeData.h"

namespace orbital_core {

	struct EarthOrientation::Impl {
		using C04 = earth_orientation_detail::C04Data;

		EarthOrientationProfile profile;
		std::optional<C04> c04;
		earth_orientation_detail::DeltaTData deltaT;

		static std::optional<C04> loadC04(const std::filesystem::path& directory, EarthOrientationProfile profile)
		{
			switch (profile) {
				case EarthOrientationProfile::IAU1980:
				case EarthOrientationProfile::FullC04:
					return std::optional<C04>(std::in_place, directory, profile);
			}
			throw std::invalid_argument("Unknown Earth orientation profile");
		}

		Impl(const std::filesystem::path& directory, EarthOrientationProfile selected,
		     const std::filesystem::path& orientationDirectory)
		    : profile(selected),
		      c04(loadC04(orientationDirectory.empty() ? directory / "earth_orientation" : orientationDirectory,
		                  selected)),
		      deltaT((orientationDirectory.empty() ? directory / "earth_orientation" : orientationDirectory) /
		             "ET-UT.dat")
		{
		}
	};

	EarthOrientation::EarthOrientation(const std::filesystem::path& directory, EarthOrientationProfile profile,
	                                   const std::filesystem::path& orientationDirectory)
	    : m_impl(std::make_unique<Impl>(directory, profile, orientationDirectory))
	{
	}

	EarthOrientation::~EarthOrientation() = default;

	EarthOrientationSample EarthOrientation::lookupMjdTt(double mjdTt) const noexcept
	{
		return m_impl->c04 ? m_impl->c04->lookupMjdTt(mjdTt) : EarthOrientationSample{};
	}

	bool EarthOrientation::coversObservedMjdTt(double mjdTt) const noexcept
	{
		return m_impl->c04 && m_impl->c04->coversObserved(mjdTt);
	}

	bool EarthOrientation::ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept
	{
		const auto* c04 = m_impl->c04 ? &*m_impl->c04 : nullptr;
		double seconds;
		if (c04 && c04->observedTtMinusUt1(mjdTt, seconds)) {
			mjdUt1 = mjdTt - seconds / 86400.0;
			return true;
		}
		return m_impl->deltaT.ttToUt1Mjd(mjdTt, mjdUt1);
	}

	bool EarthOrientation::ut1ToTtSeconds(double mjdUt1, double& seconds) const noexcept
	{
		const auto* c04 = m_impl->c04 ? &*m_impl->c04 : nullptr;
		if (c04) {
			double tt = mjdUt1 + 70.0 / 86400.0;
			if (!c04->coversObserved(tt))
				tt = mjdUt1 + 40.0 / 86400.0;
			for (int i = 0; i < 4; ++i) {
				double offset;
				if (!c04->observedTtMinusUt1(tt, offset))
					break;
				const double next = mjdUt1 + offset / 86400.0;
				if (next == tt) {
					seconds = offset;
					return true;
				}
				tt = next;
			}
			if (c04->observedTtMinusUt1(tt, seconds)) {
				return true;
			}
		}
		return m_impl->deltaT.offsetSeconds(mjdUt1, seconds);
	}

	EarthOrientationProfile EarthOrientation::profile() const noexcept { return m_impl->profile; }

	bool EarthOrientation::isAvailable() const noexcept { return m_impl->c04 && m_impl->c04->isAvailable(); }

	std::size_t EarthOrientation::size() const noexcept { return m_impl->c04 ? m_impl->c04->size() : 0; }

	const EarthOrientation& earthOrientation()
	{
		static const EarthOrientation table(getExecutableDir(), EarthOrientationProfile::IAU1980,
		                                    earthC04DataDirectory());
		return table;
	}

}  // namespace orbital_core
