#pragma once

#include "strewn_core/StrewnJob.h"

namespace strewn {

	enum class MaterialPreset { OrdinaryChondrite, CarbonaceousChondrite, Iron, Custom };

	// Only material ranges change; entry states, H, numerical controls and maps stay intact.
	// OC uses the strong branch. The table does not specify a weak-branch mixture weight.
	STREWN_CORE_EXPORT void applyMaterialPreset(JobSpec& job, MaterialPreset preset);
	STREWN_CORE_EXPORT MaterialPreset materialPreset(const JobSpec& job);

	// Accepts only a complete, single-encounter Earth covariance MC at 100 km.
	// Uses its saved entries directly; never repeats orbital sampling/propagation.
	STREWN_CORE_EXPORT JobSpec prepareAtmosphericJob(const std::filesystem::path& mcArchive);

	// Fixes the convergence grid from entry trajectories using the same coverage
	// policy as automatic weather/terrain preparation. Statistical settings survive.
	STREWN_CORE_EXPORT ConvergenceConfig atmosphericConvergenceGrid(std::span<const EntryState> entries,
	                                                                ConvergenceConfig options = {});

	// Round-trippable configuration for file/fixed-state jobs, including all defaults.
	// EntrySamples must first be archived by their producer.
	STREWN_CORE_EXPORT void saveJobConfiguration(const std::filesystem::path& path, const JobSpec& job);
	STREWN_CORE_EXPORT void validate(const MonteCarloConfig& options);

}  // namespace strewn
