#pragma once

#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <string>
#include <string_view>
#include <vector>

namespace orbital_core::obs {

	class ORBITAL_CORE_API Vfcc17Rules {
	   public:
		[[nodiscard]] bool load(std::string_view text, std::string& error);
		[[nodiscard]] std::size_t size() const noexcept;
		[[nodiscard]] ObservationUncertainty lookup(std::string_view station, char technique, char catalog,
		                                            char program, double jdUtc) const noexcept;

	   private:
		struct Rule {
			std::string station, techniques, catalogs, programs;
			double startJd = -1.0e10;
			double endJd = 1.0e10;
			double sigmaRa = 1.0;
			double sigmaDec = 1.0;
		};

		std::vector<Rule> m_rules;
	};

}  // namespace orbital_core::obs
