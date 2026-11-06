#ifndef OBS_NODE_H
#define OBS_NODE_H

#include "orbital_core/obs/ParsedObservation.h"
#include "orbital_core/ObservationIdentity.h"
#include <string>
#include <optional>
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"

#include <iostream>

struct ORBITAL_CORE_API ObservatoryRecord {
	std::string code;
	double longitude = 0.0;
	double cosPhi = 0.0;
	double sinPhi = 0.0;
	std::string name;
};

class ORBITAL_CORE_API ObsNode {
   public:
	ObsNode();
	ObsNode(const orbital_core::obs::ParsedObservation& data, const ObservatoryRecord& observatory);
	~ObsNode() = default;

	void print(std::ostream& os = std::cout) const noexcept;

	orbital_core::ObservationId id{orbital_core::newObservationId()};
	std::optional<orbital_core::obs::RadarObservation> radar;

	[[nodiscard]] int measurementDimensions() const noexcept { return radar ? radar->dimensions() : 2; }

	double JD;
	// JD and its UTC-named split fields retain the submitted time scale.
	orbital_core::time::Scale timeScale = orbital_core::time::Scale::UTC;
	// Throws for unavailable historical clock data, rather than fitting a
	// fabricated epoch. Callers retain their existing split-MJD arithmetic.
	[[nodiscard]] double ttOffsetSeconds() const;
	// Dynamical reception epoch, converted from the declared source scale.
	// Split input timestamps retain their MJD precision; JD is only a view.
	[[nodiscard]] double epochMjdTdb() const;

	[[nodiscard]] double epochJdTdb() const { return orbital_core::time::jdFromMjd(epochMjdTdb()); }

	int utcMjdDay = 0;
	double utcSecondsOfDay = 0.0;
	bool hasSplitUtc = false;
	double RA;
	double DEC;
	double MAG;
	double longitude;
	double cosPhi;
	double sinPhi;
	double LST;
	char note1 = ' ';
	char note2 = ' ';
	char programCode = ' ';
	std::string adesProgramCode;
	std::string notes;
	bool isInterstellar = false;
	char magBand = ' ';  // MPC photometric band, col 71 in 80-column records
	std::optional<bool> nuclearMagnitude;
	std::string magString;  // Submitted numeric magnitude text, without the band.
	double rmsMag = 0.0;    // Observer-supplied photometric sigma (ADES).
	std::optional<orbital_core::obs::ImportedObservationFit> importedFit;
	char catalogCode = ' ';
	std::string astCat;
	orbital_core::obs::AstrometricConvention astrometricConvention{
	    orbital_core::obs::AstrometricConvention::CatalogReferenced};
	// Observer-supplied astrometric uncertainty from ADES (rmsRA / rmsDec /
	// rmsCorr), arcsec.  Absent from 80-column records, and absent from ADES
	// records whose subFmt says they originated as 80-column ones, so the
	// weighting table remains the fallback. Used when useObserverSigmas is on.
	bool hasObserverSigma = false;
	double rmsRaArcsec = 0.0;
	double rmsDecArcsec = 0.0;
	double rmsCorr = 0.0;
	// Observer timing uncertainty, seconds.  rmsTime is the random component
	// the fit weights with; uncTime is the semi-systematic one, carried but
	// not yet modelled.
	bool hasRmsTime = false;
	double rmsTimeSec = 0.0;
	bool hasUncTime = false;
	double uncTimeSec = 0.0;

	// Catalogue debias actually applied to RA/DEC above, arcsec, in the same
	// units the correction was applied in: dAlpha (already divided by cos DEC
	// inside AstrometryBias::correction) and dDelta.  RA/DEC are stored
	// CORRECTED; submittedCoordinates preserves the original angles and
	// ecliptic direction exactly, without a trigonometric round-trip on restore.
	// Zero when debiasing is off or the catalogue has
	// no correction (Gaia).  Note ESA's .rwo uses the opposite convention:
	// raw coordinates plus this column.
	double biasRaArcsec = 0.0;
	double biasDecArcsec = 0.0;
	bool hasTimeAccuracy = false;
	double timeAccuracyDays = 1.0e-6;
	double raAccuracyArcsec = 1.5e-2;
	double decAccuracyArcsec = 1.0e-2;
	std::string obsCode;
	std::string obsName;
	Vec3D obsVector;
	Vec3D vector;

	bool isSatellite = false;
	Vec3D satelliteOffsetEclAu{0.0, 0.0, 0.0};
	std::optional<orbital_core::obs::RovingObserver> rovingObserver;

	struct Coordinates {
		double raRad;
		double decRad;
		Vec3D directionEcl;
	};

	std::optional<Coordinates> submittedCoordinates;

	[[nodiscard]] double submittedRaRad() const noexcept
	{
		return submittedCoordinates ? submittedCoordinates->raRad : RA;
	}

	[[nodiscard]] double submittedDecRad() const noexcept
	{
		return submittedCoordinates ? submittedCoordinates->decRad : DEC;
	}
};

#endif  /// OBS_NODE_H
