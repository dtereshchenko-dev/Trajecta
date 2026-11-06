#include "AstrometricSigma.h"
#include "Vfcc17Rules.h"

#include "FileReader.h"
#include "logger/Logger.h"
#include "orbital_core/Constants.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"

#include <filesystem>
#include <mutex>
#include <string>

namespace orbital_core::obs {
	namespace {

		std::once_flag g_vfcc17Once;
		Vfcc17Rules g_vfcc17Rules;

		void loadVfcc17Rules()
		{
			FileReader reader;
			const auto path = getExecutableDir() / "weights" / "vfcc17.rules";
			reader.read(path.wstring().c_str());
			if (reader.error() != FileReader::Error::NoError) {
				LOG_ERROR << "[Weights] cannot read " << path.string() << ": " << reader.errorString()
				          << "; falling back to 1 arcsec";
				return;
			}
			std::string error;
			if (!g_vfcc17Rules.load(reader.data(), error)) {
				LOG_ERROR << "[Weights] " << error << "; falling back to 1 arcsec";
				return;
			}
			LOG_INFO << "[Weights] loaded " << g_vfcc17Rules.size() << " vfcc17 rules from " << path.string();
		}

	}  // namespace

	ORBITAL_CORE_API void ensureSigmaTable(AstrometricWeightProfile profile)
	{
		switch (profile) {
			case AstrometricWeightProfile::Vfcc17:
				std::call_once(g_vfcc17Once, loadVfcc17Rules);
				break;
			case AstrometricWeightProfile::Uniform:
				break;
		}
	}

	ORBITAL_CORE_API ObservationUncertainty lookupUncertainty(const ObsNode& obs, double jd,
	                                                          AstrometricWeightProfile profile)
	{
		ensureSigmaTable(profile);

		switch (profile) {
			case AstrometricWeightProfile::Vfcc17:
				return g_vfcc17Rules.lookup(obs.obsCode, obs.note2, obs.catalogCode, obs.programCode, jd);

			case AstrometricWeightProfile::Uniform:
				return ObservationUncertainty{
				    .sigmaRaRad = 0.5 * Constants::ARCSEC_RAD,
				    .sigmaDecRad = 0.5 * Constants::ARCSEC_RAD,
				    .resolved = true,
				};
		}

		LOG_ERROR << "[Weights] unsupported weight profile " << static_cast<int>(profile)
		          << "; falling back to 1 arcsec";
		return ObservationUncertainty{};
	}

	ORBITAL_CORE_API bool usesObserverSigma(const ObsNode& obs, bool enabled) noexcept
	{
		if (!obs.hasObserverSigma || !(obs.rmsRaArcsec > 0.0) || !(obs.rmsDecArcsec > 0.0) || obs.note2 == 'X' ||
		    obs.note2 == 'x')
			return false;
		return enabled || obs.note2 == 'E';
	}

}  // namespace orbital_core::obs
