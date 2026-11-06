#include "McBoundaryScan.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "orbital_core/hazard/PlanetData.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard::mc {
	void scan(McClone& clone, const IntervalMcConfig& config, NBodyPropagator& propagator, IEphemerisSource& ephemeris,
	          std::stop_token cancel)
	{
		constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
		const auto target = toEphemerisPlanet(config.target);
		std::vector<double> targets;
		if (config.startMjdTdb > clone.epochMjdTdb)
			targets.push_back(config.startMjdTdb);
		targets.push_back(config.endMjdTdb);
		std::vector<NBodyPropagator::StateVel> outputs;
		const auto& x = clone.state;
		const auto terminal = propagator.propagateUntilSurfaceMjd(
		    {x[0], x[1], x[2]}, {x[3], x[4], x[5]}, clone.epochMjdTdb, targets, outputs, target,
		    planetRadiusKm(config.target), ephemeris, clone.nonGrav,
		    config.target == Planet::Earth ? NBodyPropagator::SurfaceShape::ReferenceEllipsoid
		                                   : NBodyPropagator::SurfaceShape::Sphere,
		    config.altitudeKm, cancel);
		if (cancel.stop_requested())
			throw std::runtime_error("MC cancelled");
		if (terminal.hitSurface && terminal.mjd == clone.epochMjdTdb)
			throw std::runtime_error("MC clone starts inside the requested boundary");
		if (!terminal.hitSurface && terminal.mjd != config.endMjdTdb)
			throw std::runtime_error("Incomplete MC interval propagation");
		if (!terminal.hitSurface)
			return;
		clone.entryBeforeInterval = terminal.mjd < config.startMjdTdb;
		if (clone.entryBeforeInterval)
			return;
		const auto body = ephemeris.queryMjd(target, terminal.mjd);
		const auto r = (terminal.R - body.R) * auKm;
		const auto v = (terminal.V - body.V) * (auKm / 86400.0);
		if (!std::isfinite(r.dot(r)) || !std::isfinite(v.dot(v)))
			throw std::runtime_error("Non-finite MC boundary state");
		const auto fixed = icrfToBodyFixedMjd(config.target, r, terminal.mjd, &ephemeris);
		const auto ll = latLonFromBodyFixed(config.target, fixed);
		if (!ll.valid)
			throw std::runtime_error("MC body-fixed conversion failed");
		clone.entry = McBoundaryEntry{terminal.mjd, r, v, ll.latitudeDeg, ll.longitudeDeg};
	}

	void summarize(IntervalMcResult& result)
	{
		std::vector<double> epochs;
		for (const auto& clone : result.clones) {
			switch (clone.status) {
				case McCloneStatus::Complete: {
					++result.completed;
					result.entriesBeforeInterval += clone.entryBeforeInterval;
					if (!clone.entry)
						break;
					++result.entryClones;
					const double epoch = clone.entry->mjdTdb;
					if (!std::isfinite(epoch))
						throw std::runtime_error("Non-finite MC entry epoch");
					epochs.push_back(epoch);
					break;
				}
				case McCloneStatus::FitFailed:
					++result.failedFit;
					break;
				case McCloneStatus::PropagationFailed:
					++result.failedPropagation;
					break;
				case McCloneStatus::Cancelled:
					++result.cancelledClones;
					break;
				case McCloneStatus::Unstarted:
					++result.unstarted;
					break;
			}
		}
		std::sort(epochs.begin(), epochs.end());
		for (const double epoch : epochs) {
			if (result.events.empty() || epoch - result.events.back().endMjdTdb > eventMaxGapDays)
				result.events.push_back({.startMjdTdb = epoch, .endMjdTdb = epoch});
			auto& row = result.events.back();
			row.endMjdTdb = epoch;
			++row.entryClones;
		}
		for (auto& row : result.events) {
			const double n = result.completed, p = row.entryClones / n, z = 1.959963984540054;
			const double denominator = 1 + z * z / n;
			const double center = (p + z * z / (2 * n)) / denominator;
			const double half = z * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator;
			row.probability = p;
			row.probabilityLow = std::max(0.0, center - half);
			row.probabilityHigh = std::min(1.0, center + half);
		}
	}
}  // namespace orbital_core::hazard::mc
