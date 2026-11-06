#pragma once

#include "strewn_core/Environment.h"
#include <span>

namespace strewn {

	struct FragmentGroup {
		double massKg{};
		double multiplicity{};
	};

	struct FragmentBatch {
		std::vector<Fragment> children;
		double dustKg{};
		double residualKg{};
	};

	class STREWN_CORE_EXPORT FragmentationModel {
	   public:
		explicit FragmentationModel(PhysicsConfig physics);
		// Interrupted grouping/splitting throws; no partial child set is returned.
		std::vector<FragmentGroup> groups(
		    const Fragment& parent, const RunParameters& parameters, std::stop_token cancel = {},
		    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max()) const;
		// One [0, 1) draw per group; explicit draws permit paired numerical tests.
		FragmentBatch split(
		    const Fragment& parent, const RunParameters& parameters, const AirSample& air,
		    std::span<const FragmentGroup> groups, std::span<const double> draws, std::stop_token cancel = {},
		    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max()) const;

	   private:
		PhysicsConfig m_physics;
	};

}  // namespace strewn
