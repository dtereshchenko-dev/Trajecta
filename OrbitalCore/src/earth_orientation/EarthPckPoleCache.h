#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"

#include <filesystem>
#include <vector>

namespace orbital_core {
	// Constructed while CSPICE is serialized; lookup only reads immutable data.
	class ORBITAL_CORE_API EarthPckPoleCache {
	   public:
		bool addFile(const std::filesystem::path& path);
		bool pole(double mjdTdb, Vec3D& result) const noexcept;

	   private:
		struct Segment {
			double start = 0;
			double end = 0;
			double initial = 0;
			double interval = 0;
			std::size_t recordSize = 0;
			std::size_t recordCount = 0;
			std::vector<double> data;
		};

		std::vector<Segment> m_segments;
		bool m_valid = true;
	};
}  // namespace orbital_core
