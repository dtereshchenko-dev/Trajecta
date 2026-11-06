#pragma once

#include "obs/IObservationReader.h"
#include <string>
#include <string_view>
#include <vector>

// ---------------------------------------------------------------------------
// The legacy MPC 80-column ASCII record, plus the whitespace-tokenised variant
// people hand-assemble from MPEC circulars.
//
// Everything the format knows stays here: column offsets, the note-1 and
// note-2 letters, the single-letter catalogue code, the '1'/'2' unit flag on
// satellite continuation lines, and the reported-precision conventions.  The
// parsing itself is delegated to MPECParser, whose intermediate record mirrors
// the wire format one field per column; this reader converts that into the
// neutral ParsedObservation.
//
// The format's cost is measurable, and is the reason ADES exists: RA is
// rounded to 0.001s, DEC to 0.01" and the timestamp to 1e-6 d, and there is no
// field at all for the observer's own uncertainty.
// ---------------------------------------------------------------------------

namespace orbital_core::obs {

	class ORBITAL_CORE_API Mpc80ColumnReader final : public IObservationReader {
	   public:
		Mpc80ColumnReader() = default;

		[[nodiscard]] const char* formatName() const noexcept override;

		// Loosest sniff of the registered readers, hence last in priority: any
		// line long enough to hold the record with a declination sign where the
		// format puts it, or — for the tokenised variant — enough whitespace
		// fields to carry a date, a position and an observatory code.
		[[nodiscard]] bool canRead(std::string_view text) const noexcept override;

		[[nodiscard]] std::vector<ParsedObservation> read(const std::string& text, std::string& error) const override;
	};

}  // namespace orbital_core::obs
