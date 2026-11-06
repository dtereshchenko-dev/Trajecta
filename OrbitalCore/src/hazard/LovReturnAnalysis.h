#pragma once

#include <span>
#include <vector>

namespace orbital_core::hazard {

	struct LovTraceEntry {
		int recordIndex{-1};
		int sampleIndex{-1};
		double jd{};
	};

	struct LovReturnAnalysisConfig {
		double showerSpanDays{90.0};
		double showerGapDays{30.0};
	};

	struct LovReturnFilament {
		std::vector<int> entryIndices;
	};

	[[nodiscard]] std::vector<LovReturnFilament> buildLovReturnFilaments(std::span<const LovTraceEntry> entries,
	                                                                     const LovReturnAnalysisConfig& cfg);

}  // namespace orbital_core::hazard
