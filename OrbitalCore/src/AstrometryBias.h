#pragma once

#include "orbital_core/AstrometricWeights.h"
#include <memory>
#include <vector>
#include "orbital_core/OrbitalCoreExport.h"
#include <filesystem>

class ObsNode;

class ORBITAL_CORE_API AstrometryBias {
   public:
	struct Correction {
		double dRA_arcsec{};
		double dDEC_arcsec{};
	};

	AstrometryBias();
	~AstrometryBias();

	AstrometryBias(const AstrometryBias&) = delete;
	AstrometryBias& operator=(const AstrometryBias&) = delete;
	AstrometryBias(AstrometryBias&&) noexcept;
	AstrometryBias& operator=(AstrometryBias&&) noexcept;

	void load(const std::filesystem::path& biasDatPath);

	[[nodiscard]] bool loaded() const noexcept;
	[[nodiscard]] std::size_t rowCount() const noexcept;
	[[nodiscard]] std::size_t catalogCount() const noexcept;
	[[nodiscard]] Correction correction(double raRad, double decRad, char catalogCode, double jd) const;

	static constexpr int NSIDE = 64;

	[[nodiscard]] static int ang2pixRing(double theta, double phi, int nside = NSIDE) noexcept;

   private:
	class Impl;
	std::unique_ptr<Impl> p_impl;
};

namespace Calculations {
	enum class AstrometryBiasScheme { Efcc18, Fcct14 };
	[[nodiscard]] std::filesystem::path astrometryBiasPathForScheme(const std::filesystem::path& runtimeDirectory,
	                                                                AstrometryBiasScheme scheme);
	[[nodiscard]] ORBITAL_CORE_API std::filesystem::path astrometryBiasPathForProfile(
	    const std::filesystem::path& runtimeDirectory, orbital_core::AstrometricWeightProfile profile);

	ORBITAL_CORE_API
	std::size_t debiasObservations(std::vector<ObsNode>& obs, const std::filesystem::path& biasPath);
}  // namespace Calculations
