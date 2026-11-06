#include "strewn_core/FragmentationModel.h"
#include "strewn_core/TrajectoryIntegrator.h"
#include "VectorMath.h"
#include "RunDeadline.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <stdexcept>
#include <utility>

namespace strewn {
	namespace {

		void validateParent(const Fragment& parent)
		{
			if (!std::isfinite(parent.state.massKg) || parent.state.massKg <= 0.0 ||
			    !std::isfinite(parent.strengthPa) || parent.strengthPa <= 0.0 || !std::isfinite(parent.multiplicity) ||
			    parent.multiplicity <= 0.0 || !std::isfinite(parent.timeS) || parent.timeS < 0.0 ||
			    !detail::finite(parent.state.positionEcefM) || !detail::finite(parent.state.velocityEcefMps) ||
			    parent.generation == std::numeric_limits<std::uint32_t>::max())
				throw std::invalid_argument("Invalid fragmentation parent");
		}

	}  // namespace

	FragmentationModel::FragmentationModel(PhysicsConfig physics) : m_physics(std::move(physics))
	{
		validate(m_physics, {});
	}

	std::vector<FragmentGroup> FragmentationModel::groups(const Fragment& parent, const RunParameters& parameters,
	                                                      std::stop_token cancel,
	                                                      std::chrono::steady_clock::time_point deadline) const
	{
		validate(parameters);
		validateParent(parent);
		detail::checkInterruption(cancel, deadline);
		const double dust = m_physics.dustAtEveryBreakup || parent.generation == 0 ? parameters.dustFraction : 0.0;
		const double available = parent.state.massKg * (1.0 - dust);
		if (available < m_physics.minMassKg)
			return {};
		const std::size_t phase = parent.generation == 0 ? 0 : 1;
		const auto individual = m_physics.individualFragments[phase];
		const double bins = m_physics.binsPerDecade[phase];
		std::vector<FragmentGroup> result;
		result.reserve(std::min<std::size_t>(individual, 100) + 40);
		std::map<int, FragmentGroup> tail;
		std::size_t count = 0;
		const auto add = [&](double mass) {
			if (mass < m_physics.minMassKg)
				return;
			if (count++ < individual) {
				result.push_back({mass, 1.0});
			} else {
				auto& bin = tail[static_cast<int>(std::floor(std::log10(mass) * bins))];
				bin.massKg += mass;
				bin.multiplicity += 1.0;
			}
		};
		const double largest = parameters.largestFraction * available;
		add(largest);
		double left = available - largest;
		const double inverse =
		    parameters.beta * std::pow(parameters.largestFraction, 1.0 - parameters.beta) / (1.0 - parameters.beta);
		const double origin = std::pow(parameters.largestFraction, -parameters.beta);
		for (std::size_t k = 1; left > 0.0; ++k) {
			if ((k & 255) == 0)
				detail::checkInterruption(cancel, deadline);
			if (k == std::numeric_limits<std::size_t>::max())
				throw std::overflow_error("Fragment mass index overflow");
			const double mass = std::pow(origin + static_cast<double>(k) * inverse, -1.0 / parameters.beta) * available;
			if (mass < m_physics.minMassKg)
				break;
			if (mass >= left) {
				if (m_physics.retainMassRemainder)
					add(left);
				break;
			}
			add(mass);
			left -= mass;
		}
		for (auto [index, bin] : tail) {
			detail::checkInterruption(cancel, deadline);
			bin.massKg /= bin.multiplicity;
			const auto representatives =
			    std::min(m_physics.representativesPerBin, static_cast<std::size_t>(bin.multiplicity));
			const double weight = bin.multiplicity / static_cast<double>(representatives);
			for (std::size_t i = 0; i < representatives; ++i) {
				if ((i & 255) == 0)
					detail::checkInterruption(cancel, deadline);
				result.push_back({bin.massKg, i + 1 == representatives
				                                  ? bin.multiplicity - weight * (representatives - 1)
				                                  : weight});
			}
		}
		detail::checkInterruption(cancel, deadline);
		return result;
	}

	FragmentBatch FragmentationModel::split(const Fragment& parent, const RunParameters& parameters,
	                                        const AirSample& air, std::span<const FragmentGroup> grouped,
	                                        std::span<const double> draws, std::stop_token cancel,
	                                        std::chrono::steady_clock::time_point deadline) const
	{
		validate(parameters);
		validateParent(parent);
		detail::checkInterruption(cancel, deadline);
		if (!std::isfinite(air.densityKgM3) || air.densityKgM3 <= 0.0 || !detail::finite(air.windEnuMps))
			throw std::invalid_argument("Invalid air at breakup");
		if (draws.size() != grouped.size())
			throw std::invalid_argument("Fragment kick draw count differs from group count");
		FragmentBatch result;
		const double dust = m_physics.dustAtEveryBreakup || parent.generation == 0 ? parameters.dustFraction : 0.0;
		result.dustKg = parent.state.massKg * parent.multiplicity * dust;
		result.residualKg = parent.state.massKg * parent.multiplicity - result.dustKg;
		if (grouped.empty())
			return result;
		const auto position = geodetic(parent.state.positionEcefM);
		const auto wind = enuToEcef(air.windEnuMps, position.latitudeRad, position.longitudeRad);
		Vector3 relative, unit;
		const double speed = detail::norm(parent.state.velocityEcefMps);
		if (speed == 0.0)
			throw std::invalid_argument("Cannot fragment a stationary parent");
		for (std::size_t i = 0; i < 3; ++i) {
			relative[i] = parent.state.velocityEcefMps[i] - parameters.windScale * wind[i];
			unit[i] = parent.state.velocityEcefMps[i] / speed;
		}
		auto e1 = detail::cross(unit, std::abs(unit[2]) < 0.99 ? Vector3{0, 0, 1} : Vector3{1, 0, 0});
		const double e1Length = detail::norm(e1);
		for (double& value : e1)
			value /= e1Length;
		const auto e2 = detail::cross(unit, e1);
		const double parentDiameter = std::cbrt(6.0 * parent.state.massKg / (detail::pi * parameters.densityKgM3));
		const double relativeSpeed = detail::norm(relative);
		double totalMass = 0.0;
		Vector3 mean{};
		result.children.reserve(grouped.size());
		for (std::size_t i = 0; i < grouped.size(); ++i) {
			if ((i & 255) == 0)
				detail::checkInterruption(cancel, deadline);
			const auto& group = grouped[i];
			if (!std::isfinite(draws[i]) || draws[i] < 0.0 || draws[i] >= 1.0 || !std::isfinite(group.massKg) ||
			    group.massKg <= 0.0 || !std::isfinite(group.multiplicity) || group.multiplicity <= 0.0)
				throw std::invalid_argument("Invalid fragment group or kick draw");
			const double diameter = std::cbrt(6.0 * group.massKg / (detail::pi * parameters.densityKgM3));
			const double ratio = parentDiameter / diameter;
			const double magnitude = std::sqrt(m_physics.lateralCoefficient * 1.5 * air.densityKgM3 /
			                                   parameters.densityKgM3 * ratio * ratio) *
			                         relativeSpeed;
			const double phi = 2.0 * detail::pi * draws[i];
			Fragment child;
			child.state.positionEcefM = parent.state.positionEcefM;
			child.state.massKg = group.massKg;
			child.timeS = parent.timeS;
			child.strengthPa =
			    std::min(m_physics.strengthCapPa,
			             parent.strengthPa * std::pow(parent.state.massKg / group.massKg, m_physics.weibullAlpha));
			child.multiplicity = parent.multiplicity * group.multiplicity;
			child.parentId = parent.id;
			child.generation = parent.generation + 1;
			child.shape = m_physics.fragmentShape;
			const double mass = group.massKg * group.multiplicity;
			totalMass += mass;
			for (std::size_t c = 0; c < 3; ++c) {
				child.state.velocityEcefMps[c] = magnitude * (std::cos(phi) * e1[c] + std::sin(phi) * e2[c]);
				mean[c] += mass * child.state.velocityEcefMps[c];
			}
			result.children.push_back(child);
		}
		if (!std::isfinite(totalMass) || totalMass > parent.state.massKg * (1.0 - dust) * (1.0 + 1e-12))
			throw std::invalid_argument("Fragment groups exceed the available parent mass");
		for (std::size_t i = 0; i < result.children.size(); ++i) {
			if ((i & 255) == 0)
				detail::checkInterruption(cancel, deadline);
			auto& child = result.children[i];
			for (std::size_t c = 0; c < 3; ++c)
				child.state.velocityEcefMps[c] += parent.state.velocityEcefMps[c] - mean[c] / totalMass;
		}
		result.residualKg = std::max(0.0, result.residualKg - totalMass * parent.multiplicity);
		detail::checkInterruption(cancel, deadline);
		return result;
	}

}  // namespace strewn
