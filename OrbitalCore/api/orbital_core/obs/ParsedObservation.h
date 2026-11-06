#pragma once
#include "orbital_core/obs/RovingObserver.h"
#include "orbital_core/obs/AstrometricConvention.h"

#include <optional>
#include <string>
#include "orbital_core/TimeScale.h"
#include "orbital_core/obs/RadarObservation.h"

// ---------------------------------------------------------------------------
// Format-neutral optical or radar observation, produced by readers and consumed by
// ObsNode.  It deliberately does NOT mirror any wire format:
// Radar channels are optional native-unit measurements; optical coordinates and
// photometric fields are not applicable when radar is present.
//
//   * coordinates are radians, stored once — no sexagesimal text, no second
//     representation with a "which one is valid" flag;
//   * the observation technique is an enum, not the MPC note-2 letter;
//   * observer-offset units are an enum, not the '1'/'2' flag from the MPC
//     satellite continuation line.
//
// The predecessor (`MPECInitData`) was an 80-column record in struct form, so
// the ADES reader had to transliterate itself into MPC column conventions
// before the rest of the pipeline would accept it.  Anything a reader must
// know about a foreign format now stays inside that reader.
//
// Two fields are deliberately still textual:
//
//   catalogCode  the single-letter MPC star-catalogue code.  This is not a
//                format detail — it is the vocabulary the debias tables and
//                the weighting rules are keyed on, so it is the canonical
//                identifier here as well.  ADES carries a name ("Gaia2") and
//                its reader maps it; the 80-column reader reads the letter
//                directly.
//
//   magnitudeText  retained verbatim rather than parsed to a double, because
//                the .rwo export reproduces the reported field as printed and
//                re-formatting from a double changes trailing digits.  ObsNode
//                derives the numeric value from it.
// ---------------------------------------------------------------------------

namespace orbital_core::obs {
	enum class Technique : int {
		Unknown = 0,
		Ccd,
		Cmos,
		Photographic,
		Video,
		Encoder,
		Micrometer,
		Meridian,
		Occultation,
		NormalPlace,
		Photometer,
		CorrectedCcd,
		PhotographicB1950,
		ReplacedCcd,
		ReplacedPhotographic,
	};

	enum class OffsetUnits : int {
		None = 0,
		Au,
		Km,
	};

	// Technique -> the single MPC note-2 letter.  Like catalogCode this is shared
	// vocabulary rather than a format detail: the astrometric weighting rules and
	// the .rwo export are keyed on the letter, so it has to be reachable from the
	// neutral record.  Codes that describe observer geometry rather than technique
	// ('S' satellite, 'V' roving, 'R' radar) are not produced here — geometry is
	// carried by isSatellite, rovingObserver and radar.
	[[nodiscard]] inline char mpcNote2(Technique technique) noexcept
	{
		switch (technique) {
			case Technique::Ccd:
				return 'C';
			case Technique::Cmos:
				return 'B';
			case Technique::Photographic:
				return 'P';
			case Technique::Video:
				return 'n';
			case Technique::Encoder:
				return 'e';
			case Technique::Micrometer:
				return 'M';
			case Technique::Meridian:
				return 'T';
			case Technique::Occultation:
				return 'E';
			case Technique::NormalPlace:
				return 'N';
			case Technique::Photometer:
				return 'H';
			case Technique::CorrectedCcd:
				return 'c';
			case Technique::PhotographicB1950:
				return 'A';
			case Technique::ReplacedCcd:
				return 'X';
			case Technique::ReplacedPhotographic:
				return 'x';
			case Technique::Unknown:
				break;
		}
		return ' ';
	}

	// Techniques whose note-2 letter is kept for a space-based observer position.
	// ADES occultation astrometry is referred to a geocentric point (sys/pos), yet
	// it is weighted and timed as an occultation ('E') despite its space-based
	// geometry; a deprecated record stays deprecated ('X').
	[[nodiscard]] inline bool techniqueOverridesGeometry(Technique technique) noexcept
	{
		return technique == Technique::Occultation || technique == Technique::ReplacedCcd ||
		       technique == Technique::ReplacedPhotographic;
	}

	// Archived RWO diagnostics, not manual include/exclude instructions for a new
	// fit. Optical bias/residuals are sky-plane arcseconds; radar uses channel units.
	struct ImportedObservationFit {
		bool selected{true};
		double bias1{}, bias2{};
		bool forceWeight1{}, forceWeight2{};
		std::optional<double> residual1, residual2;
	};

	struct ParsedObservation {
		std::string permanentDesignation;
		std::string classificationDesignation;
		std::optional<ImportedObservationFit> importedFit;
		std::string id;
		std::optional<RadarObservation> radar;
		// ---- submitted time ----
		// `jd` is always populated.  The split form mirrors reference and keeps tens
		// of microseconds that would be lost by adding a fractional day to a full
		// Julian date; readers that have it set `hasSplitUtc`. The legacy UTC-named
		// split fields use timeScale too: MPC/RWO before 1972 are UT1, including
		// legacy MPC records wrapped in ADES; native ADES is UTC.
		double jd = 0.0;
		time::Scale timeScale = time::Scale::UTC;
		int utcMjdDay = 0;
		double utcSecondsOfDay = 0.0;
		bool hasSplitUtc = false;

		// ---- astrometry, radians ----
		double raRad = 0.0;
		double decRad = 0.0;

		// ---- identity ----
		std::string objectName;
		std::string obsCode;

		// ---- classification ----
		char catalogCode = ' ';
		std::string astCat;
		AstrometricConvention astrometricConvention{AstrometricConvention::CatalogReferenced};
		std::optional<double> gaiaPhotocenterRadiusKm;
		Technique technique = Technique::Unknown;
		char programCode = ' ';       // Single-character MPC program used by weight rules.
		std::string adesProgramCode;  // Original base-62 prog, including codes beyond MPC80.
		std::string notes;            // ADES observing notes, not a program identifier.
		bool isInterstellar = false;

		// ---- photometry ----
		std::string magnitudeText;
		char magBand = ' ';
		std::optional<bool> nuclearMagnitude;
		double rmsMag = 0.0;  // ADES photometric uncertainty; zero when absent.

		// ---- precision as reported by the submitter ----
		// raAccuracyArcsec is on the sky, i.e. it already carries cos(DEC).
		// Granularity of the submitted timestamp.  The 80-column record always
		// states it by field width; ADES restates it in precTime.  `hasTimeAccuracy`
		// separates "the submitter said 1e-6 d" from the struct default.
		bool hasTimeAccuracy = false;
		double timeAccuracyDays = 1.0e-6;
		double raAccuracyArcsec = 1.5e-2;
		double decAccuracyArcsec = 1.0e-2;

		// ---- uncertainty supplied by the observer (ADES only) ----
		bool hasObserverSigma = false;
		double rmsRaArcsec = 0.0;
		double rmsDecArcsec = 0.0;
		double rmsCorr = 0.0;
		// rmsTime is the observer's RANDOM timing uncertainty; uncTime is the
		// semi-systematic one. Both are in seconds. Only rmsTime enters the
		// independent measurement covariance; uncTime is retained for export
		// fidelity and for a future correlated-timing model.
		bool hasRmsTime = false;
		double rmsTimeSec = 0.0;
		bool hasUncTime = false;
		double uncTimeSec = 0.0;

		// ---- per-observation terrestrial position ----
		std::optional<RovingObserver> rovingObserver;

		// ---- space-based observer ----
		// Offset from the geocentre to the observer, in the J2000 equatorial
		// frame.  ObsNode converts to ecliptic AU.
		bool isSatellite = false;
		OffsetUnits offsetUnits = OffsetUnits::None;
		double offsetX = 0.0;
		double offsetY = 0.0;
		double offsetZ = 0.0;

		// ---- provenance, for logging and for the .rwo header ----
		std::string sourceFormat;

		// Duplicate detection during IOD pool assembly: two records describe the
		// same measurement when the epoch, the object, the position and the
		// station agree.  Provenance and reported precision are deliberately not
		// compared — the same observation re-served in another format is still
		// the same observation.
		[[nodiscard]] bool operator==(const ParsedObservation& other) const noexcept
		{
			return jd == other.jd && radar == other.radar && objectName == other.objectName && raRad == other.raRad &&
			       decRad == other.decRad && magnitudeText == other.magnitudeText && obsCode == other.obsCode;
		}
	};

}  // namespace orbital_core::obs
