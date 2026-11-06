#pragma once

#include "orbital_core/Constants.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/obs/RadarObservation.h"
#include "orbital_core/obs/AstrometricConvention.h"

#include <string>

namespace orbital_core::obs {
	// Per-observatory geodetic data parsed from `OBSCODE.dat` (longitude in
	// degrees, cosPhi/sinPhi = parallax constants).  Spacecraft entries
	// (HST=250, WISE=C51, TESS=C57, …) carry zero geodetic fields — for
	// satellites the actual observer position lives in
	// `PreparedObservation::observerOffsetEclAu`, populated at prepare()
	// time from `ObsNode::satelliteOffsetEclAu` (the raw 's' continuation).
	struct ObserverData {
		double longitude{};
		double cosPhi{};
		double sinPhi{};
		std::string obsCode;
	};

	// Pre-resolved observation tuple consumed by the LM hot path:
	// optical fields are precomputed so optical residual evaluation costs no
	// ephemeris lookup, σ-table search, GAST/EOP query, or nutation
	// series evaluation. Consumed by OrbitFitModel and arc extension.
	//
	//   jdTdb/mjdTdb - reception epoch in TDB for both optical and radar models.
	//   jdUtc  — Raw UTC observation epoch (logging / MPC export only).
	//   raObs, decObs — observed RA/DEC, radians.
	//   earthR, earthV — Earth heliocentric state at jdTdb (DE440-equivalent).
	//   observer — geodetic data for ground stations (zero for satellites).
	//   observerOffsetEclAu — geocenter→observer vector in J2000 ecliptic
	//                          frame, AU.  Unified for ground + satellites:
	//                            - satellites: directly from MPC 's' line
	//                              (already in ecliptic AU at parse time).
	//                            - ground stations: precomputed in
	//                              ObservationPipeline::prepare() via the
	//                              full GAST → TIRS → polar-motion → CIRS →
	//                              precession-nutation → J2000 ecliptic
	//                              chain, depending only on (jdTdb, observer)
	//                              which are fixed for an observation.
	//                          `toTopocentric` uses `R_helio - earthR - observerOffsetEclAu`
	//                          with catalogue-referenced solar light deflection
	//                          (LightDeflection.h) for catalogue-referenced optics.
	//   sigmaRa, sigmaDec — per-obs astrometric σ in radians (1 arcsec default).
	//                       The optimizer decides whether to apply whitening.
	//   catalogCode — MPC col-72 catalog letter ('*' wildcard match key).
	struct PreparedObservation {
		double jdTdb{};
		// Authoritative reception epoch without JD rounding.
		double mjdTdb{};
		double jdUtc{};
		double mjdUtc{};
		double raObs{};
		double decObs{};
		// Tangent-plane metric of the submitted coordinates, before debiasing.
		double raMetricCosDec{1.0};
		Vec3D earthR;
		Vec3D earthV;
		ObserverData observer;
		std::optional<RadarObservation> radar;
		ObserverData transmitter;

		[[nodiscard]] int measurementDimensions() const noexcept { return radar ? radar->dimensions() : 2; }

		bool isSatellite{false};
		Vec3D observerOffsetEclAu{0.0, 0.0, 0.0};
		Vec3D observerVelocityEclAuPerDay{0.0, 0.0, 0.0};

		double sigmaRa{Constants::ARCSEC_RAD};
		double sigmaDec{Constants::ARCSEC_RAD};
		bool forceWeightRa{false};
		bool forceWeightDec{false};
		// Observer-supplied ADES rmsRA/rmsDec replaced the table sigma. Such
		// observations still count toward the vfcc17 same-station batch size
		// but are not relaxed by it. This is distinct from the forceWeight*
		// flags exported in RWO; observer sigmas alone do not force those flags.
		bool observerSigma{false};
		// Sigma-table timing uncertainty in days; zero means absent.
		double timeSigmaDays{0.0};
		// RA/DEC correlation reported by the observer (ADES rmsCorr), in
		// [-1, 1].  Zero leaves the error ellipse axis-aligned.
		double sigmaCorr{0.0};
		double mag{0.0};
		char magBand{' '};
		std::optional<bool> nuclearMagnitude;
		double magSigma{0.0};
		char note1{' '};
		char note2{' '};
		char catalogCode{' '};
		AstrometricConvention astrometricConvention{AstrometricConvention::CatalogReferenced};
		double photocenterRadiusAu{0.0};
	};
}  // namespace orbital_core::obs
