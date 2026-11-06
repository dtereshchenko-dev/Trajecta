#pragma once

#include "strewn_core/Environment.h"

namespace strewn {

	STREWN_CORE_EXPORT void validate(const PhysicsConfig& config, const IntegrationOptions& integration);
	STREWN_CORE_EXPORT void validate(const RunParameters& parameters);
	STREWN_CORE_EXPORT double projectedArea(double massKg, double densityKgM3, Shape shape);
	STREWN_CORE_EXPORT double dragCoefficient(double mach);

	class STREWN_CORE_EXPORT TrajectoryIntegrator {
	   public:
		TrajectoryIntegrator(const Environment& environment, const PhysicsConfig& physics,
		                     const IntegrationOptions& options = {});
		TrajectoryIntegrator(const Environment&&, const PhysicsConfig&, const IntegrationOptions& = {}) = delete;
		FlightResult advance(const Fragment& fragment, const RunParameters& parameters,
		                     std::stop_token cancel = {}) const;
		// Diagnostic observations of accepted states; does not change the step sequence.
		FlightResult advance(
		    const Fragment& fragment, const RunParameters& parameters, std::stop_token cancel,
		    const FlightObserver& observer,
		    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max()) const;

	   private:
		const Environment& m_environment;
		PhysicsConfig m_physics;
		IntegrationOptions m_options;
	};

}  // namespace strewn
