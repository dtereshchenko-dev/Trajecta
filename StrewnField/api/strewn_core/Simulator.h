#pragma once

#include "strewn_core/TrajectoryIntegrator.h"
#include <span>

namespace strewn {

	class STREWN_CORE_EXPORT Simulator {
	   public:
		// Environment is immutable and must outlive this simulator and its calls.
		Simulator(const Environment& environment, PhysicsConfig physics = {}, IntegrationOptions integration = {});
		Simulator(const Environment&&, PhysicsConfig = {}, IntegrationOptions = {}) = delete;
		RunResult runOne(const EntryState& entry, const RunParameters& parameters, std::uint64_t seed,
		                 std::uint64_t runId = 0, std::stop_token cancel = {}, double runTimeoutS = 300.0) const;
		// Replays one realization with diagnostic traces, without sampling a new ensemble.
		RunResult runOne(const EntryState& entry, const RunParameters& parameters, std::uint64_t seed,
		                 std::uint64_t runId, std::stop_token cancel, const FragmentObserver& observer,
		                 double runTimeoutS = 300.0) const;
		// Callbacks are serialized on the calling thread, in run-ID order.
		RunSummary run(std::span<const EntryState> entries, const MonteCarloConfig& options,
		               const RunControl& control = {}, const RunConsumer& consume = {}) const;

	   private:
		const Environment& m_environment;
		PhysicsConfig m_physics;
		IntegrationOptions m_integration;
	};

}  // namespace strewn
