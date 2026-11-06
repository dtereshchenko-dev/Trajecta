#pragma once

#include "logger/Logger.h"
#include "EnvironmentVariable.h"

#include <cstdlib>
#include <iomanip>
#include <locale>
#include <sstream>

// Per-force acceleration audit, the propagation-side counterpart of
// hazard::audit. Off unless ORBITAL_FORCE_AUDIT=1, so the hot path keeps a
// single predictable branch and normal runs are byte-identical.
//
//   ORBITAL_FORCE_AUDIT=1              enable
//   ORBITAL_FORCE_AUDIT_COUNT=n        dump the first n force evaluations of each
//                                      armed phase (default 1)
//   ORBITAL_FORCE_AUDIT_ENCOUNTER=n    also dump n evaluations inside the Earth
//                                      J2 region (default 0); the Earth
//                                      orientation terms are identically zero
//                                      outside it, so first-step records alone
//                                      say nothing about pole/nutation/polar motion.
namespace orbital_core::forces::audit {

	inline int setting(const char* name, int fallback) noexcept
	{
		const auto valueText = environmentVariable(name);
		if (!valueText || valueText->empty())
			return fallback;
		const char* raw = valueText->c_str();
		char* end = nullptr;
		const long value = std::strtol(raw, &end, 10);
		if (end == raw)
			return fallback;
		return static_cast<int>(value);
	}

	inline bool enabled() noexcept
	{
		static const bool value = setting("ORBITAL_FORCE_AUDIT", 0) == 1;
		return value;
	}

	inline long& openingBudget() noexcept
	{
		static thread_local long value = enabled() ? setting("ORBITAL_FORCE_AUDIT_COUNT", 1) : 0;
		return value;
	}

	inline long& encounterBudget() noexcept
	{
		static thread_local long value = enabled() ? setting("ORBITAL_FORCE_AUDIT_ENCOUNTER", 0) : 0;
		return value;
	}

	// Which run phase the records belong to. Budgets are granted per phase, so
	// the differential correction cannot spend what the survey needs: the fit
	// and the LOV survey both open at the same epoch and would otherwise be
	// indistinguishable in the log.
	inline const char*& phase() noexcept
	{
		static thread_local const char* value = "fit";
		return value;
	}

	// Re-grants both budgets and relabels subsequent records. Reached from
	// outside the library only through the exported arm() in
	// api/orbital_core/forces/ForceAuditControl.h: these thread_local budgets
	// live in the DLL, and a header-inline copy linked into the executable
	// would be a different object entirely.
	inline void armInternal(const char* name) noexcept
	{
		if (!enabled())
			return;
		phase() = name;
		openingBudget() = setting("ORBITAL_FORCE_AUDIT_COUNT", 1);
		encounterBudget() = setting("ORBITAL_FORCE_AUDIT_ENCOUNTER", 0);
	}

	// Sequence number of the force evaluation being dumped, so that records
	// from the opening steps and from an encounter never get interleaved.
	inline long& evaluationIndex() noexcept
	{
		static thread_local long value = 0;
		return value;
	}

	// Claims one record. `inEncounter` lets encounter evaluations draw on
	// their own budget once the opening budget is spent.
	inline bool claim(bool inEncounter) noexcept
	{
		if (!enabled())
			return false;
		if (openingBudget() > 0) {
			--openingBudget();
			return true;
		}
		if (inEncounter && encounterBudget() > 0) {
			--encounterBudget();
			return true;
		}
		return false;
	}

	template <class Writer>
	void write(const char* stage, Writer&& writer) noexcept
	{
		try {
			std::ostringstream out;
			out.imbue(std::locale::classic());
			out << std::setprecision(17) << "[FORCE_AUDIT] stage=" << stage << " phase=" << phase()
			    << " eval=" << evaluationIndex();
			writer(out);
			LOG_INFO << out.str();
		} catch (...) {
			// Diagnostics must never disturb the propagation they measure.
		}
	}
}  // namespace orbital_core::forces::audit
