#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/CloseApproach.h"

#include <Eigen/Core>

#include <span>
#include <vector>

namespace orbital_core::hazard {

	struct ORBITAL_CORE_API BPlaneBasis {
		Vec3D s{};
		Vec3D xi{};
		Vec3D zeta{};
		bool valid{false};
	};

	struct ORBITAL_CORE_API BPlanePoint {
		double xi{};
		double zeta{};
	};

	struct ORBITAL_CORE_API BPlaneCovariance {
		double meanXi{};
		double meanZeta{};
		double sigmaXi{};
		double sigmaZeta{};
		double correlation{};
		int numSamples{};
		bool analytic{false};

		[[nodiscard]] bool isValid() const noexcept
		{
			return (analytic || numSamples > 1) && sigmaXi > 0.0 && sigmaZeta > 0.0;
		}
	};

	struct ORBITAL_CORE_API BPlaneJacobian {
		Eigen::Vector2d bScaledKm{Eigen::Vector2d::Zero()};
		Eigen::Matrix<double, 2, 6> dBdXTcaKm{Eigen::Matrix<double, 2, 6>::Zero()};
		BPlaneBasis basis{};
		double gamma{1.0};
		bool valid{false};
	};

	struct ORBITAL_CORE_API ScaledBPlanePoint {
		Eigen::Vector2d bScaledKm{Eigen::Vector2d::Zero()};
		BPlaneBasis basis{};
		double gamma{1.0};
		double vInfKmS{0.0};
		bool valid{false};
	};

	[[nodiscard]] ORBITAL_CORE_API BPlaneBasis computeBPlaneBasis(const CloseApproachEvent& event) noexcept;

	// Build a b-plane basis directly from a relative-velocity vector
	// (clone-velocity − planet-velocity at TCA, AU/day).  Used by
	// observation-noise MC: each clone has its OWN v_∞ direction, so
	// projecting clone position onto the NOMINAL basis gives a
	// mixed-frame number — for short-arc encounters the per-clone
	// basis is required to recover the actual b/γ impact metric.
	[[nodiscard]] ORBITAL_CORE_API BPlaneBasis computeBPlaneBasisFromVrel(const Vec3D& vRelAuPerDay) noexcept;

	[[nodiscard]] ORBITAL_CORE_API BPlanePoint mapToBPlane(const Vec3D& missVecKm, const BPlaneBasis& basis) noexcept;

	[[nodiscard]] ORBITAL_CORE_API BPlaneCovariance fitBPlaneCovariance(std::span<const Vec3D> samplesKm,
	                                                                    const BPlaneBasis& basis) noexcept;

	[[nodiscard]] ORBITAL_CORE_API double gravitationalFocusing(double vInfKmS, Planet planet) noexcept;

	[[nodiscard]] ORBITAL_CORE_API ScaledBPlanePoint computeScaledBPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS,
	                                                                          Planet planet) noexcept;
}  // namespace orbital_core::hazard
