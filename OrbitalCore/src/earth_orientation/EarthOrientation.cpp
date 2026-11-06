#include "EarthOrientation.h"
#include <filesystem>
#include "EopData.h"
#include "orbital_core/Types.h"

namespace orbital_core {

	struct EarthOrientation::Impl {
		using C04 = earth_orientation_detail::C04Data;

		EarthOrientationProfile profile;
		std::optional<C04> c04;
		earth_orientation_detail::FinalsData finals;
		earth_orientation_detail::DeltaTData deltaT;

		static std::optional<C04> loadC04(const std::filesystem::path& directory, EarthOrientationProfile profile)
		{
			switch (profile) {
				case EarthOrientationProfile::IAU1980:
				case EarthOrientationProfile::FullC04:
					return std::optional<C04>(std::in_place, directory / "earth_orientation", profile);
				case EarthOrientationProfile::IAU2000:
					return std::nullopt;
			}
			throw std::invalid_argument("Unknown Earth orientation profile");
		}

		Impl(const std::filesystem::path& directory, EarthOrientationProfile selected)
		    : profile(selected),
		      c04(loadC04(directory, selected)),
		      finals(directory / "earth_orientation" / "finals2000A.all"),
		      deltaT(directory / "earth_orientation" / "ET-UT.dat")
		{
		}
	};

	EarthOrientation::EarthOrientation(const std::filesystem::path& directory, EarthOrientationProfile profile)
	    : m_impl(std::make_unique<Impl>(directory, profile))
	{
	}

	EarthOrientation::~EarthOrientation() = default;

	EarthOrientationSample EarthOrientation::lookupMjdTt(double mjdTt) const noexcept
	{
		return m_impl->c04 ? m_impl->c04->lookupMjdTt(mjdTt) : m_impl->finals.lookupMjdTt(mjdTt);
	}

	bool EarthOrientation::ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept
	{
		return m_impl->finals.ttToUt1Mjd(mjdTt, mjdUt1) || m_impl->deltaT.ttToUt1Mjd(mjdTt, mjdUt1);
	}

	bool EarthOrientation::ut1ToTtSeconds(double mjdUt1, double& seconds) const noexcept
	{
		return m_impl->finals.ut1ToTtSeconds(mjdUt1, seconds) || m_impl->deltaT.offsetSeconds(mjdUt1, seconds);
	}

	EarthOrientationProfile EarthOrientation::profile() const noexcept { return m_impl->profile; }

	bool EarthOrientation::isAvailable() const noexcept
	{
		return m_impl->c04 ? m_impl->c04->isAvailable() : m_impl->finals.isAvailable();
	}

	std::size_t EarthOrientation::size() const noexcept
	{
		return m_impl->c04 ? m_impl->c04->size() : m_impl->finals.size();
	}

	const EarthOrientation& earthOrientation()
	{
		static const EarthOrientation table(getExecutableDir());
		return table;
	}

}  // namespace orbital_core
