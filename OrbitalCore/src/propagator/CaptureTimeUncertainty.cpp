#include "orbital_core/RelativeMotion.h"

#include "RelativeMotionInternal.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"
#include "uncertainty/LinearizedCovariance.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace orbital_core {
	RelativeMotionTimeUncertainty relativeMotionTimeUncertainty(const RelativeMotionRequest& request,
	                                                            const std::vector<RelativeMotionSample>& samples,
	                                                            const std::vector<TemporaryCaptureInterval>& captures,
	                                                            const std::vector<HillSpherePassage>& hillPassages,
	                                                            const FitRunControl& control)
	{
		RelativeMotionTimeUncertainty result{std::vector<CaptureTimeUncertainty>(captures.size()),
		                                     std::vector<CaptureTimeUncertainty>(hillPassages.size())};
		const auto& solution = request.solution;
		const auto factor =
		    uncertainty::covarianceFactor(solution.covariance, 6 + solution.nonGrav.fitParameterCount());
		if (!solution.forces.enabled || samples.size() < 2 || (captures.empty() && hillPassages.empty()) || !factor)
			return result;

		struct Boundary {
			std::size_t event, left;
			bool start, energy, hill;
			double fraction, slope;
		};

		std::vector<Boundary> boundaries;
		std::vector<std::size_t> indices;
		const auto addBoundary = [&](std::size_t event, bool start, double epoch) {
			const auto right = std::lower_bound(samples.begin(), samples.end(), epoch,
			                                    [](const auto& sample, double value) { return sample.mjdTdb < value; });
			if (right == samples.begin() || right == samples.end())
				return;
			const std::size_t left = static_cast<std::size_t>(right - samples.begin() - 1);
			const auto& a = samples[left];
			const auto& b = samples[left + 1];
			const double dt = b.mjdTdb - a.mjdTdb;
			if (!(dt > 0))
				return;
			const double values[2][2] = {{a.specificEnergyAu2PerDay2, b.specificEnergyAu2PerDay2},
			                             {a.distanceAu - 3 * a.hillRadiusAu, b.distanceAu - 3 * b.hillRadiusAu}};
			int selected = -1;
			double best = std::numeric_limits<double>::infinity(), fraction = 0, slope = 0;
			for (int k = 0; k < 2; ++k) {
				const double x = values[k][0], y = values[k][1];
				if (!std::isfinite(x) || !std::isfinite(y) || x == y || (x > 0 && y > 0) || (x < 0 && y < 0))
					continue;
				const double f = x / (x - y);
				const double error = std::abs((a.mjdTdb + f * dt) - epoch);
				if (selected >= 0 && std::abs(error - best) < dt * 1e-8)
					return;
				if (error < best) {
					selected = k;
					best = error;
					fraction = f;
					slope = (y - x) / dt;
				}
			}
			if (selected < 0 || best > std::max(1e-9, dt * 1e-7) || !std::isfinite(slope) || slope == 0)
				return;
			boundaries.push_back({event, left, start, selected == 0, false, fraction, slope});
			indices.push_back(left);
			indices.push_back(left + 1);
		};
		for (std::size_t i = 0; i < captures.size(); ++i) {
			if (!captures[i].startClipped)
				addBoundary(i, true, captures[i].startMjdTdb);
			if (!captures[i].endClipped)
				addBoundary(i, false, captures[i].endMjdTdb);
		}
		const auto addHillBoundary = [&](std::size_t event, bool entry, double epoch) {
			const auto right = std::lower_bound(samples.begin(), samples.end(), epoch,
			                                    [](const auto& sample, double value) { return sample.mjdTdb < value; });
			const std::size_t index = static_cast<std::size_t>(right - samples.begin());
			for (const std::size_t left : {index > 0 ? index - 1 : index, index}) {
				if (left + 1 >= samples.size())
					continue;
				const auto& a = samples[left];
				const auto& b = samples[left + 1];
				const double dt = b.mjdTdb - a.mjdTdb;
				const double x = a.distanceAu - a.hillRadiusAu;
				const double y = b.distanceAu - b.hillRadiusAu;
				if (!(dt > 0) || !std::isfinite(x) || !std::isfinite(y) || !(entry ? x >= 0 && y < 0 : x <= 0 && y > 0))
					continue;
				const double fraction = x / (x - y);
				const double slope = (y - x) / dt;
				if (!std::isfinite(slope) || slope == 0 ||
				    std::abs(a.mjdTdb + fraction * dt - epoch) > std::max(1e-9, dt * 1e-7))
					continue;
				boundaries.push_back({event, left, entry, false, true, fraction, slope});
				indices.push_back(left);
				indices.push_back(left + 1);
				return;
			}
		};
		for (std::size_t i = 0; i < hillPassages.size(); ++i) {
			if (!hillPassages[i].entryClipped)
				addHillBoundary(i, true, hillPassages[i].entryMjdTdb);
			if (!hillPassages[i].exitClipped)
				addHillBoundary(i, false, hillPassages[i].exitMjdTdb);
		}
		if (boundaries.empty() || control.cancel.stop_requested())
			return result;
		std::sort(indices.begin(), indices.end());
		indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
		std::vector<double> epochs;
		epochs.reserve(indices.size());
		for (auto i : indices)
			epochs.push_back(samples[i].mjdTdb);
		auto* source = getDefaultEphemerisSource();
		if (!source || !source->isAvailable())
			return result;
		auto propagator = detail::relativeMotionPropagator(solution.forces, *source);
		const auto states =
		    propagator.propagateAllWithStmAndSensMjd(solution.position, solution.velocity, solution.epochMjdTdb, epochs,
		                                             solution.nonGrav, nullptr, control.cancel);
		if (control.cancel.stop_requested() || states.size() != indices.size())
			return result;
		for (std::size_t i = 0; i < indices.size(); ++i) {
			const auto& sample = samples[indices[i]];
			const auto planet = source->queryMjd(request.planet, sample.mjdTdb);
			const double difference = (states[i].R - planet.R - sample.positionAu).magnitude();
			if (!std::isfinite(difference) || difference > std::max(1e-8, sample.distanceAu * 1e-5))
				return result;
		}
		double mu = 0;
		for (const auto& body : forces::detail::kPerturbers)
			if (body.planet == request.planet)
				mu = body.gm;
		if (!(mu > 0))
			return result;
		const auto gradient = [&](std::size_t sampleIndex, bool energy) -> Eigen::RowVectorXd {
			const auto pos = static_cast<std::size_t>(std::lower_bound(indices.begin(), indices.end(), sampleIndex) -
			                                          indices.begin());
			const auto& sample = samples[sampleIndex];
			const auto& r = sample.positionAu;
			const auto& v = sample.velocityAuPerDay;
			const double radius = sample.distanceAu;
			Eigen::Matrix<double, 1, 6> g;
			const double radial = energy ? mu / (radius * radius * radius) : 1 / radius;
			g << radial * r.x, radial * r.y, radial * r.z, energy ? v.x : 0, energy ? v.y : 0, energy ? v.z : 0;
			return g * uncertainty::stateJacobian(states[pos], solution.nonGrav);
		};
		for (const auto& boundary : boundaries) {
			if (control.cancel.stop_requested())
				return result;
			const Eigen::RowVectorXd derivative = ((1 - boundary.fraction) * gradient(boundary.left, boundary.energy) +
			                                       boundary.fraction * gradient(boundary.left + 1, boundary.energy)) /
			                                      boundary.slope;
			const double sigma = 3 * (derivative * *factor).norm();
			if (std::isfinite(sigma)) {
				auto& event = (boundary.hill ? result.hillPassages : result.captures)[boundary.event];
				(boundary.start ? event.start3SigmaDays : event.end3SigmaDays) = sigma;
			}
		}
		return result;
	}

	std::vector<CaptureTimeUncertainty> captureTimeUncertainty(const RelativeMotionRequest& request,
	                                                           const std::vector<RelativeMotionSample>& samples,
	                                                           const std::vector<TemporaryCaptureInterval>& captures,
	                                                           const FitRunControl& control)
	{
		return relativeMotionTimeUncertainty(request, samples, captures, {}, control).captures;
	}
}  // namespace orbital_core
