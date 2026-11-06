#ifndef ADESPARSER_H
#define ADESPARSER_H

#include <string>
#include <string_view>
#include <vector>

#include "obs/IObservationReader.h"

// ---------------------------------------------------------------------------
// ADES — IAU Astrometry Data Exchange Standard, adopted by IAU Commission 20
// in 2015 as the successor to the 80-column ASCII record.  Reference
// implementation and the normative field list: github.com/IAU-ADES/ADES-Master
// (`xml/adesmaster.xml`); submission rules: minorplanetcenter.net/iau/info/ADES
//
// Two interchangeable serialisations carry identical content:
//
//   XML   <ades version="2022"><optical><stn>F52</stn>...</optical>...</ades>
//   PSV   pipe-separated, a keyword record naming the columns then data rows,
//         optionally preceded by '#'-prefixed header lines.
//
// Why it earns its place next to the 80-column reader: that record rounds RA
// to 0.001s, DEC to 0.01" and the timestamp to 1e-6 d, and has no field for
// the observer's own uncertainty.  Measured on the 2020 VV arc, the rounding
// alone is ~21 mas in RA, ~14 mas in DEC and ~8 ms in time, and reading ADES
// instead moved the fitted orbit 264 km closer to JPL's.
//
// This reader extracts optical, occultation and radar measurements consumed by
// the fit. subFmt=PRE/M47/M92 retains the original MPC time convention (UT1
// before 1972, UTC afterwards): conversion to ADES ISO text does not change
// the time scale. Native ADES timestamps remain UTC. Other submission metadata
// and telescope descriptions not needed by the fit are ignored.
//
// Occultations (<occultation>, or PSV rows carrying raStar/decStar) use
// Dec = decStar + deltaDec,
// RA = raStar + deltaRA / cos(Dec), the reference point from sys/pos1..3 as
// for a space-based observer, technique 'E'.  The stated precision of
// deltaRA/deltaDec is the accuracy, and replaces rmsRA/rmsDec when it exceeds
// both.  dist/pa offsets are rejected rather than converted.  deprecated=X
// marks a replaced observation ('X'), as the 80-column record does; vfcc17
// gives it the 999" exclusion weight.
// Unknown fields are skipped rather than rejected, so a
// batch carrying extra columns still loads.
// ---------------------------------------------------------------------------

namespace orbital_core::ades {

	class ORBITAL_CORE_API AdesReader final : public obs::IObservationReader {
	   public:
		AdesReader() = default;

		[[nodiscard]] const char* formatName() const noexcept override;
		[[nodiscard]] bool canRead(std::string_view text) const noexcept override;
		[[nodiscard]] std::vector<obs::ParsedObservation> read(const std::string& text,
		                                                       std::string& error) const override;
	};

}  // namespace orbital_core::ades

#endif  // ADESPARSER_H
