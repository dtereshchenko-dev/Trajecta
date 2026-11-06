#ifndef RESIDUALS_H
#define RESIDUALS_H

#include "orbital_core/OrbitalCoreExport.h"

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

struct ResidualWeight {
	// Marginal uncertainties after weighting/timing; RA is RA*cos(DEC).
	double sigmaRaArcsec{1.0};
	double sigmaDecArcsec{1.0};
	// Measurement-covariance norm, before leverage correction.
	double chiSquared{0.0};
	bool forceWeightRa{false};
	bool forceWeightDec{false};
	bool weighted{false};
	// Uses the formal orbit covariance, with minus/plus for fitted/excluded rows.
	std::optional<double> postFitChiSquared;
};

struct Residual {
	struct RadarChannel {
		double observed{};
		double computed{};
		double residual{};
		double sigma{};
	};

	struct Radar {
		std::optional<RadarChannel> delaySeconds;
		std::optional<RadarChannel> dopplerHz;
		std::string transmitter;
		double frequencyHz{};
		double chiSquared{};
		std::optional<double> postFitChiSquared;
	};

	double jd;
	double raObs;         // observed RA (radians)
	double decObs;        // observed DEC (radians)
	double dRA_arcsec;    // O-C in RA*cos(DEC), arcseconds
	double dDEC_arcsec;   // O-C in DEC, arcseconds
	double total_arcsec;  // sqrt(dRA^2 + dDEC^2)
	std::string obsCode;
	double mag{0.0};
	double helioDistAu{0.0};
	double geoDistAu{0.0};
	double phaseAngleRad{0.0};
	char magBand{' '};
	std::optional<bool> nuclearMagnitude;
	double magSigma{0.0};
	std::optional<ResidualWeight> astrometricWeight;
	std::optional<Radar> radar;
};

class ORBITAL_CORE_API Residuals {
   public:
	Residuals() = default;
	~Residuals() = default;

	void add(const Residual& r);
	void clear();

	double rms() const;
	size_t count() const noexcept;
	const std::vector<Residual>& data() const noexcept;

	void print() const;
	void print(std::ostream& os) const;
	bool save(const std::filesystem::path& filePath) const;

   private:
	void printHeader(std::ostream& os) const;
	std::vector<Residual> m_residuals;
};

#endif  // RESIDUALS_H
