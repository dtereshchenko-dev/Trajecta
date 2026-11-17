#pragma once

#include "NBodyPropagator.h"
#include "TwoBodyPropagation.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <stop_token>
#include <utility>
#include <vector>

namespace orbital_core::obs {

	struct LightTimeSeries {
		std::vector<double> emissionMjd;
		std::vector<double> emissionOffsetDays;
		std::vector<NBodyPropagator::StateVel> states;
	};

	struct DenseLightTimeSeries {
		LightTimeSeries light;
		std::vector<IntegratedStateFull> variationalStates;
	};

	inline Vec3D sunPositionMjd(double mjd)
	{
		if (auto* source = getDefaultEphemerisSource())
			if (const auto state = source->barycentricSunMjd(mjd))
				return state->R;
		return {};
	}

	inline Vec3D sunVelocityMjd(double mjd)
	{
		if (auto* source = getDefaultEphemerisSource())
			if (const auto state = source->barycentricSunMjd(mjd))
				return state->V;
		return {};
	}

	inline Vec3D sunAccelerationMjd(double mjd)
	{
		if (auto* source = getDefaultEphemerisSource()) {
			constexpr double step = 1e-3;
			const auto before = source->barycentricSunMjd(mjd - step);
			const auto after = source->barycentricSunMjd(mjd + step);
			if (before && after)
				return (after->V - before->V) / (2.0 * step);
		}
		return {};
	}

	inline Vec3D apparentEmissionPosition(const NBodyPropagator::StateVel& state, double emissionMjd,
	                                      double receptionMjd, double emissionOffsetDays)
	{
		return state.R + sunPositionMjd(emissionMjd) + sunVelocityMjd(emissionMjd) * emissionOffsetDays -
		       sunPositionMjd(receptionMjd);
	}

	template <typename Sample>
	DenseLightTimeSeries solveDenseLightTimeSeries(Sample&& sample, const std::vector<double>& receptionMjd,
	                                               const std::vector<Vec3D>& observers,
	                                               const std::vector<bool>& optical)
	{
		const std::size_t count = receptionMjd.size();
		if (observers.size() != count || optical.size() != count)
			throw std::invalid_argument("Light-time input size mismatch");
		DenseLightTimeSeries result;
		result.light.emissionMjd = receptionMjd;
		result.light.emissionOffsetDays.resize(count);
		result.light.states.resize(count);
		result.variationalStates.resize(count);
		for (std::size_t i = 0; i < count; ++i) {
			const auto reception = sample(receptionMjd[i]);
			if (!reception)
				throw std::runtime_error("Reception epoch outside dense trajectory");
			if (!optical[i]) {
				result.variationalStates[i] = *reception;
				result.light.states[i] = {reception->R, reception->V};
				continue;
			}
			const Vec3D sun = sunPositionMjd(receptionMjd[i]);
			double epoch = receptionMjd[i] - (reception->R - observers[i]).magnitude() / Constants::AU_PER_DAY;
			bool converged = false;
			for (int iteration = 0; iteration < 8; ++iteration) {
				const auto state = sample(epoch);
				if (!state)
					throw std::runtime_error("Emission epoch outside dense trajectory");
				const double distance = (state->R + sunPositionMjd(epoch) - sun - observers[i]).magnitude();
				if (!std::isfinite(distance) || !(distance > 0.0))
					throw std::runtime_error("Invalid light-time distance");
				const double next = receptionMjd[i] - distance / Constants::AU_PER_DAY;
				const double ulp = std::nextafter(epoch, std::numeric_limits<double>::infinity()) - epoch;
				if (std::abs(next - epoch) <= ulp) {
					const double offset = (receptionMjd[i] - epoch) - distance / Constants::AU_PER_DAY;
					result.light.emissionMjd[i] = epoch;
					result.light.emissionOffsetDays[i] = offset;
					result.light.states[i] = {state->R + state->V * offset, state->V};
					result.variationalStates[i] = *state;
					converged = true;
					break;
				}
				epoch = next;
			}
			if (!converged)
				throw std::runtime_error("Light-time iteration did not converge");
		}
		return result;
	}

	inline DenseLightTimeSeries propagateDenseLightTimeSeries(const Vec3D& initialR, const Vec3D& initialV,
	                                                          double epochMjd, const std::vector<double>& epochs,
	                                                          const std::vector<Vec3D>& observers,
	                                                          const std::vector<bool>& optical,
	                                                          const NBodyPropagator& propagator, const NonGravModel& ng,
	                                                          bool variational, bool sensitivities)
	{
		if (observers.size() != epochs.size() || optical.size() != epochs.size())
			throw std::invalid_argument("Light-time input size mismatch");
		if (epochs.empty())
			return {};
		const auto [first, last] = std::minmax_element(epochs.begin(), epochs.end());
		double observerRadius = 0.0;
		for (const auto& observer : observers)
			observerRadius = std::max(observerRadius, observer.magnitude());
		const double span = std::max(std::abs(*first - epochMjd), std::abs(*last - epochMjd));
		const double margin =
		    std::max(0.25, 2.0 * (initialR.magnitude() + initialV.magnitude() * span + observerRadius + 1.0) /
		                       Constants::AU_PER_DAY);
		double start = *first - margin;
		if (auto* source = getDefaultEphemerisSource()) {
			if (source->barycentricSunMjd(*first) && !source->barycentricSunMjd(start))
				start = *first;
		}
		auto trajectory =
		    propagator.propagateDenseMjd(initialR, initialV, epochMjd, start, *last, ng, variational, sensitivities);
		auto sample = [&](double mjd) {
			auto state = trajectory.at(mjd);
			if (!state && mjd < start) {
				start = mjd - 0.25;
				propagator.extendDenseBackwardMjd(trajectory, start, ng, variational, sensitivities);
				state = trajectory.at(mjd);
			}
			return state;
		};
		return solveDenseLightTimeSeries(sample, epochs, observers, optical);
	}

	// Solve all emission epochs together so one integration covers the complete arc.
	inline LightTimeSeries solveLightTimeSeries(const Vec3D& initialR, const Vec3D& initialV, double epochMjd,
	                                            const std::vector<double>& receptionMjd,
	                                            const std::vector<NBodyPropagator::StateVel>& receptionStates,
	                                            const std::vector<Vec3D>& observerPositions,
	                                            const std::vector<bool>& optical, bool nBody,
	                                            const NBodyPropagator& propagator, const NonGravModel& ng,
	                                            std::stop_token cancel = {})
	{
		const std::size_t count = receptionMjd.size();
		if (receptionStates.size() != count || observerPositions.size() != count || optical.size() != count)
			throw std::invalid_argument("Light-time input size mismatch");
		LightTimeSeries result;
		result.emissionMjd = receptionMjd;
		result.emissionOffsetDays.resize(count);
		std::vector<Vec3D> receptionSun(count);
		for (std::size_t i = 0; i < count; ++i) {
			if (!optical[i])
				continue;
			receptionSun[i] = sunPositionMjd(receptionMjd[i]);
			result.emissionMjd[i] -= (receptionStates[i].R - observerPositions[i]).magnitude() / Constants::AU_PER_DAY;
		}
		const StateVectors initial(initialR, initialV / Constants::GAUSS_K);
		for (int iteration = 0; iteration < 8; ++iteration) {
			if (cancel.stop_requested())
				return {};
			if (nBody)
				result.states = propagator.propagateAllWithVelocityMjd(initialR, initialV, epochMjd, result.emissionMjd,
				                                                       ng, cancel);
			else {
				result.states.resize(count);
				for (std::size_t i = 0; i < count; ++i) {
					const auto state = Calculations::propagateKeplerState(initial, result.emissionMjd[i] - epochMjd);
					if (!state.isValid())
						throw std::runtime_error("Light-time propagation failed");
					result.states[i] = {state.first, state.second * Constants::GAUSS_K};
				}
			}
			if (cancel.stop_requested())
				return {};
			if (result.states.size() != count)
				throw std::runtime_error("Incomplete light-time propagation");
			bool converged = true;
			std::vector<double> next = result.emissionMjd;
			for (std::size_t i = 0; i < count; ++i) {
				if (!optical[i])
					continue;
				const Vec3D ray =
				    result.states[i].R + sunPositionMjd(result.emissionMjd[i]) - receptionSun[i] - observerPositions[i];
				const double distance = ray.magnitude();
				if (!std::isfinite(distance) || !(distance > 0.0))
					throw std::runtime_error("Invalid light-time distance");
				next[i] = receptionMjd[i] - distance / Constants::AU_PER_DAY;
				result.emissionOffsetDays[i] =
				    (receptionMjd[i] - result.emissionMjd[i]) - distance / Constants::AU_PER_DAY;
				const double ulp = std::nextafter(result.emissionMjd[i], std::numeric_limits<double>::infinity()) -
				                   result.emissionMjd[i];
				if (std::abs(next[i] - result.emissionMjd[i]) > ulp)
					converged = false;
			}
			if (converged) {
				for (std::size_t i = 0; i < count; ++i)
					result.states[i].R += result.states[i].V * result.emissionOffsetDays[i];
				return result;
			}
			result.emissionMjd = std::move(next);
		}
		throw std::runtime_error("Light-time iteration did not converge");
	}

	inline Eigen::Matrix3d emissionPositionFactor(const Vec3D& ray, const Vec3D& barycentricVelocity)
	{
		const Eigen::Vector3d n(ray.x / ray.magnitude(), ray.y / ray.magnitude(), ray.z / ray.magnitude());
		const Eigen::Vector3d v(barycentricVelocity.x, barycentricVelocity.y, barycentricVelocity.z);
		return Eigen::Matrix3d::Identity() - v * n.transpose() / (Constants::AU_PER_DAY + n.dot(v));
	}
}  // namespace orbital_core::obs
