#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "orbital_core/obs/ParsedObservation.h"
#include "orbital_core/OrbitalCoreExport.h"

// ---------------------------------------------------------------------------
// One reader per observation wire format.  Each implementation owns everything
// its format knows — column offsets, catalogue naming, technique letters, unit
// flags — and hands back the format-neutral ParsedObservation.  Nothing about
// one format leaks into another or into the pipeline behind them.
//
// Adding a format means adding a reader and registering it; no call site
// changes and no enum grows.
// ---------------------------------------------------------------------------

namespace orbital_core::obs {

	class ORBITAL_CORE_API IObservationReader {
	   public:
		virtual ~IObservationReader() = default;

		IObservationReader(const IObservationReader&) = delete;
		IObservationReader& operator=(const IObservationReader&) = delete;

		// Name for logs and for the provenance line of exported products.
		[[nodiscard]] virtual const char* formatName() const noexcept = 0;

		// Content sniff over (a prefix of) the payload.  Must be cheap and must
		// not throw: every registered reader is asked in turn until one accepts.
		// Deciding by file extension is deliberately not an option — an ADES batch
		// arrives as .txt as often as not.
		[[nodiscard]] virtual bool canRead(std::string_view text) const noexcept = 0;

		// Parses the whole payload.  On failure returns an empty vector and sets
		// `error`; on success `error` is left empty.
		[[nodiscard]] virtual std::vector<ParsedObservation> read(const std::string& text,
		                                                          std::string& error) const = 0;

	   protected:
		IObservationReader() = default;
	};

	// Readers in priority order.  The MPC 80-column reader is deliberately last:
	// its sniff is the loosest, so the stricter formats get first refusal.
	[[nodiscard]] ORBITAL_CORE_API std::vector<const IObservationReader*> observationReaders() noexcept;

	// First reader whose canRead() accepts the payload, or nullptr when none does.
	[[nodiscard]] ORBITAL_CORE_API const IObservationReader* selectObservationReader(std::string_view text) noexcept;

}  // namespace orbital_core::obs
