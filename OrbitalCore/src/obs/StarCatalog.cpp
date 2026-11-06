#include "orbital_core/obs/StarCatalog.h"

#include <algorithm>

namespace orbital_core::obs {
	namespace {

		constexpr StarCatalog catalogs[]{{' ', "UNK", "Unknown"},
		                                 {'a', "USNOA1", "USNO-A1"},
		                                 {'b', "USNOSA1", "USNO-SA1"},
		                                 {'c', "USNOA2", "USNO-A2"},
		                                 {'d', "USNOSA2", "USNO-SA2"},
		                                 {'e', "UCAC1", "UCAC1"},
		                                 {'f', "Tyc1", "Tycho-1"},
		                                 {'g', "Tyc2", "Tycho-2"},
		                                 {'h', "GSC1.0", "GSC 1.0"},
		                                 {'i', "GSC1.1", "GSC 1.1"},
		                                 {'j', "GSC1.2", "GSC 1.2"},
		                                 {'k', "GSC2.2", "GSC 2.2"},
		                                 {'l', "ACT", "ACT"},
		                                 {'m', "GSCACT", "GSC-ACT"},
		                                 {'n', "SDSS8", "SDSS DR8"},
		                                 {'o', "USNOB1", "USNO-B1"},
		                                 {'p', "PPM", "PPM"},
		                                 {'q', "UCAC4", "UCAC4"},
		                                 {'r', "UCAC2", "UCAC2"},
		                                 {'s', "USNOB2", "USNO-B2"},
		                                 {'t', "PPMXL", "PPMXL"},
		                                 {'u', "UCAC3", "UCAC3"},
		                                 {'v', "NOMAD", "NOMAD"},
		                                 {'w', "CMC14", "CMC14"},
		                                 {'x', "Hip2", "Hipparcos-2"},
		                                 {'y', "Hip1", "Hipparcos-1"},
		                                 {'z', "GSC", "GSC"},
		                                 {'A', "AC", "AC"},
		                                 {'B', "SAO1984", "SAO 1984"},
		                                 {'C', "SAO", "SAO"},
		                                 {'D', "AGK3", "AGK3"},
		                                 {'E', "FK4", "FK4"},
		                                 {'F', "ACRS", "ACRS"},
		                                 {'G', "LickGas", "Lick Gaspra"},
		                                 {'H', "Ida93", "Ida 1993"},
		                                 {'I', "Perth70", "Perth 70"},
		                                 {'J', "COSMOS", "COSMOS"},
		                                 {'K', "Yale", "Yale"},
		                                 {'L', "2MASS", "2MASS"},
		                                 {'M', "GSC2.3", "GSC 2.3"},
		                                 {'N', "SDSS7", "SDSS DR7"},
		                                 {'O', "SSTRC1", "SSTRC1"},
		                                 {'P', "MPOSC3", "MPOSC3"},
		                                 {'Q', "CMC15", "CMC15"},
		                                 {'R', "SSTRC4", "SSTRC4"},
		                                 {'S', "URAT1", "URAT1"},
		                                 {'T', "URAT2", "URAT2"},
		                                 {'U', "Gaia1", "Gaia DR1"},
		                                 {'V', "Gaia2", "Gaia DR2"},
		                                 {'W', "Gaia3", "Gaia DR3"},
		                                 {'X', "Gaia3E", "Gaia EDR3"},
		                                 {'Y', "UCAC5", "UCAC5"},
		                                 {'Z', "ATLAS2", "ATLAS-2"},
		                                 {'0', "IHW", "IHW"},
		                                 {'1', "PS1_DR1", "Pan-STARRS1 DR1"},
		                                 {'2', "PS1_DR2", "Pan-STARRS1 DR2"},
		                                 {'3', "Gaia_Int", "Gaia INT"},
		                                 {'4', "GZ", "GZ"},
		                                 {'5', "UBSC", "UBSC"},
		                                 {'6', "Gaia_2016", "Gaia 2016"},
		                                 {'7', "ZZCAT", "ZZCAT"},
		                                 {'8', "APASS", "APASS"},
		                                 {'9', "UBSC", "UBSC"},
		                                 {'!', "AKARI", "AKARI"},
		                                 {'@', "AG", "AG"},
		                                 {'#', "WISE", "WISE"},
		                                 {'$', "LSST2502", "LSST2502"},
		                                 {'^', "IRASPSC", "IRASPSC"},
		                                 {'&', "SDSS9", "SDSS DR9"},
		                                 {'*', "SDSS16", "SDSS DR16"}};

		constexpr char asciiLower(char value) noexcept
		{
			return value >= 'A' && value <= 'Z' ? static_cast<char>(value + ('a' - 'A')) : value;
		}

	}  // namespace

	const StarCatalog* findStarCatalog(char mpcCode) noexcept
	{
		for (const auto& catalog : catalogs)
			if (catalog.mpcCode == mpcCode)
				return &catalog;
		return nullptr;
	}

	const StarCatalog* findStarCatalog(std::string_view adesName) noexcept
	{
		for (const auto& catalog : catalogs)
			if (std::ranges::equal(catalog.adesName, adesName, {}, asciiLower, asciiLower))
				return &catalog;
		return nullptr;
	}

}  // namespace orbital_core::obs
