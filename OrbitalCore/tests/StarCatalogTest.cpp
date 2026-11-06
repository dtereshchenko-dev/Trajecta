#include <gtest/gtest.h>
#include "orbital_core/obs/StarCatalog.h"

#include <algorithm>
#include <string>
#include <string_view>

using orbital_core::obs::findStarCatalog;

TEST(StarCatalog, CoversMpcCodesAndPreservesExistingAdesMappings)
{
	constexpr std::string_view codes = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$^&*";
	constexpr std::string_view names[]{
	    "UNK",    "USNOA1", "USNOSA1", "USNOA2",  "USNOSA2", "UCAC1",   "Tyc1",     "Tyc2",    "GSC1.0",  "GSC1.1",
	    "GSC1.2", "GSC2.2", "ACT",     "GSCACT",  "SDSS8",   "USNOB1",  "PPM",      "UCAC4",   "UCAC2",   "USNOB2",
	    "PPMXL",  "UCAC3",  "NOMAD",   "CMC14",   "Hip2",    "Hip1",    "GSC",      "AC",      "SAO1984", "SAO",
	    "AGK3",   "FK4",    "ACRS",    "LickGas", "Ida93",   "Perth70", "COSMOS",   "Yale",    "2MASS",   "GSC2.3",
	    "SDSS7",  "SSTRC1", "MPOSC3",  "CMC15",   "SSTRC4",  "URAT1",   "URAT2",    "Gaia1",   "Gaia2",   "Gaia3",
	    "Gaia3E", "UCAC5",  "ATLAS2",  "IHW",     "PS1_DR1", "PS1_DR2", "Gaia_Int", "GZ",      "UBSC",    "Gaia_2016",
	    "ZZCAT",  "APASS",  "UBSC",    "AKARI",   "AG",      "WISE",    "LSST2502", "IRASPSC", "SDSS9",   "SDSS16"};
	static_assert(codes.size() == std::size(names));
	for (std::size_t i = 0; i < codes.size(); ++i) {
		SCOPED_TRACE(codes[i]);
		const auto* catalog = findStarCatalog(codes[i]);
		ASSERT_NE(catalog, nullptr);
		EXPECT_EQ(catalog->mpcCode, codes[i]);
		EXPECT_EQ(catalog->adesName, names[i]);
		EXPECT_FALSE(catalog->name.empty());
		std::string lower(names[i]);
		std::ranges::transform(lower, lower.begin(),
		                       [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c + ('a' - 'A')) : c; });
		const auto* byName = findStarCatalog(names[i]);
		ASSERT_NE(byName, nullptr);
		EXPECT_EQ(byName, findStarCatalog(lower));
		EXPECT_EQ(byName->mpcCode, codes[i] == '9' ? '5' : codes[i]);
	}
	EXPECT_EQ(findStarCatalog('W')->name, "Gaia DR3");
	EXPECT_EQ(findStarCatalog('X')->name, "Gaia EDR3");
	EXPECT_EQ(findStarCatalog('&')->name, "SDSS DR9");
}

TEST(StarCatalog, DistinguishesUnknownNamesFromExplicitUnk)
{
	EXPECT_EQ(findStarCatalog(""), nullptr);
	EXPECT_EQ(findStarCatalog("FutureCatalog_DR9"), nullptr);
	EXPECT_EQ(findStarCatalog('?'), nullptr);
	EXPECT_EQ(findStarCatalog('\0'), nullptr);
	ASSERT_NE(findStarCatalog("uNk"), nullptr);
	EXPECT_EQ(findStarCatalog("uNk")->mpcCode, ' ');
	EXPECT_EQ(findStarCatalog('a')->adesName, "USNOA1");
	EXPECT_EQ(findStarCatalog('A')->adesName, "AC");
}
