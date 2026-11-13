#pragma once

namespace orbital_core {
	enum class RelativityMode { Full, SunOuter, Disabled };

	[[nodiscard]] constexpr const char* relativityModeName(RelativityMode mode) noexcept
	{
		switch (mode) {
			case RelativityMode::Full: return "full_eih";
			case RelativityMode::SunOuter: return "sun_outer_eih";
			case RelativityMode::Disabled: return "off";
		}
		return "unknown";
	}

	struct ForceModelSettings {
		RelativityMode relativity{RelativityMode::Full};
		bool earthJ2{true};
		bool earthJ3{true};
		bool earthJ4{true};
		bool rotatingEarthPole{true};
	};
}
