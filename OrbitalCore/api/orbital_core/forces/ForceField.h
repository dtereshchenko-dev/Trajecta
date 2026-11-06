#pragma once

#include "IForceModel.h"
#include <utility>
#include <tuple>

namespace orbital_core::forces {
	template <typename... Forces>
	    requires(ForceModel<Forces> && ...)
	class ForceField {
	   public:
		constexpr ForceField() = default;

		template <typename... Args>
		explicit constexpr ForceField(Args&&... args) : m_forces{std::forward<Args>(args)...}
		{
		}

		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& dadv) const noexcept
		{
			std::apply([&ctx, &a, &dadr, &dadv](const auto&... f) noexcept { (f.evaluate(ctx, a, dadr, dadv), ...); },
			           m_forces);
		}

		// Convenience: zero outputs first, then accumulate.
		void evaluateZeroed(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& dadv) const noexcept
		{
			a = Vec3D(0.0, 0.0, 0.0);
			dadr = Eigen::Matrix3d::Zero();
			dadv = Eigen::Matrix3d::Zero();
			evaluate(ctx, a, dadr, dadv);
		}

		// Applies fn to each force in evaluation order. Diagnostics use it to
		// accumulate through a shared context and difference the running total,
		// which is the only decomposition consistent with the inter-force
		// ephemeris cache PlanetaryPerturbations fills for EIH and EarthJ2.
		template <typename Fn>
		void forEach(Fn&& fn) const
		{
			std::apply([&fn](const auto&... f) { (fn(f), ...); }, m_forces);
		}

		[[nodiscard]] static constexpr std::size_t size() noexcept { return sizeof...(Forces); }

		[[nodiscard]] static constexpr bool empty() noexcept { return sizeof...(Forces) == 0; }

		template <typename F>
		[[nodiscard]] F& get() noexcept
		{
			return std::get<F>(m_forces);
		}

		template <typename F>
		[[nodiscard]] const F& get() const noexcept
		{
			return std::get<F>(m_forces);
		}

	   private:
		std::tuple<Forces...> m_forces;
	};
}  // namespace orbital_core::forces
