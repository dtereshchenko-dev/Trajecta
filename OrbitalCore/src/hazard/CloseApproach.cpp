#include "orbital_core/hazard/CloseApproach.h"

#include "orbital_core/Constants.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/hazard/Moid.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "NBodyPropagator.h"
#include "hazard/CloseApproachKinematics.h"
#include "hazard/HazardAudit.h"
#include "CloseApproachScan.h"
#include <stdexcept>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace orbital_core::hazard {
	namespace {
		struct Sample {
			Vec3D rTarget{};
			Vec3D vTarget{};
			Vec3D rPlanet{};
			Vec3D vPlanet{};
			double d2{};
		};

		struct SampleGrid {
			std::vector<double> grid;
			std::vector<Sample> samples;
			std::vector<NBodyPropagator::PropagationWithSurfaceStop> surfaceEntries;
		};

		[[nodiscard]] double localMoidAu(const Sample& sample, double jd) noexcept
		{
			const double invGauss = 1.0 / Constants::GAUSS_K;
			const auto targetOrbit =
			    make_orbit_from_states_vectors(StateVectors(sample.rTarget, sample.vTarget * invGauss), jd);
			const auto planetOrbit =
			    make_orbit_from_states_vectors(StateVectors(sample.rPlanet, sample.vPlanet * invGauss), jd);
			if (!targetOrbit || !planetOrbit)
				return std::numeric_limits<double>::quiet_NaN();

			const MoidResult moid = computeMoid(*targetOrbit, *planetOrbit);
			return moid.isValid() ? moid.moidAu : std::numeric_limits<double>::quiet_NaN();
		}

		void updateLocalMoid(CloseApproachEvent& event, const Sample& encounterEntry, double encounterEntryJd) noexcept
		{
			event.localMoidAu = localMoidAu(encounterEntry, encounterEntryJd);
		}

		[[nodiscard]] std::optional<SampleGrid> buildSampleGridFromState(
		    const Vec3D& r0Au, const Vec3D& v0AuPerDay, double epochTdb, const NBodyPropagator& propagator,
		    const NonGravModel* nonGrav, const CloseApproachConfig& cfg, IEphemerisSource* eph,
		    bool includeEndpoint = false, const FitRunControl* control = nullptr);

		[[nodiscard]] std::optional<SampleGrid> buildSampleGrid(const Orbit& orbit, const NBodyPropagator& propagator,
		                                                        const CloseApproachConfig& cfg, IEphemerisSource* eph)
		{
			return buildSampleGridFromState(orbit.getR(), orbit.getV() * Constants::GAUSS_K, orbit.getEpoch(),
			                                propagator, nullptr, cfg, eph);
		}

		[[nodiscard]] std::optional<SampleGrid> buildSampleGridFromState(
		    const Vec3D& r0Au, const Vec3D& v0AuPerDay, double epochTdb, const NBodyPropagator& propagator,
		    const NonGravModel* nonGrav, const CloseApproachConfig& cfg, IEphemerisSource* eph, bool includeEndpoint,
		    const FitRunControl* control)
		{
			if (cfg.horizonYears <= 0.0 || cfg.lookbackYears < 0.0 || cfg.thresholdAu <= 0.0 ||
			    cfg.encounterRadiusAu <= 0.0 || cfg.sampleStepDay <= 0.0)
				return std::nullopt;

			if (!eph)
				return std::nullopt;

			std::vector<double> grid;
			const EphemerisPlanet ephPlanet = toEphemerisPlanet(cfg.planet);

			const double t0 = epochTdb - cfg.lookbackYears * 365.25;
			const double tEnd = epochTdb + cfg.horizonYears * 365.25;
			const double dt = cfg.sampleStepDay;

			const std::size_t nSteps = static_cast<std::size_t>((tEnd - t0) / dt) + 1;

			grid.clear();
			grid.reserve(nSteps + (includeEndpoint ? 1 : 0));
			for (std::size_t i = 0; i < nSteps; ++i) {
				const double t = t0 + static_cast<double>(i) * dt;
				if (t > tEnd)
					break;

				grid.push_back(t);
			}

			if (includeEndpoint && (grid.empty() || grid.back() < tEnd))
				grid.push_back(tEnd);

			SampleGrid result;
			try {
				auto appendSample = [&](double jd, const Vec3D& r, const Vec3D& v) {
					const auto p = eph->queryJD(ephPlanet, jd);
					const Vec3D dr = r - p.R;
					result.grid.push_back(jd);
					result.samples.push_back({r, v, p.R, p.V, dr.dot(dr)});
				};
				auto runDirection = [&](const std::vector<double>& targets) {
					if (targets.empty())
						return;
					std::vector<NBodyPropagator::StateVel> states;
					std::vector<double> targetMjds;
					if (control) {
						for (double jd : targets)
							targetMjds.push_back(jd - 2400000.5);
					}
					const auto terminal =
					    control ? propagator.propagateUntilSurfaceMjd(
					                  r0Au, v0AuPerDay, epochTdb - 2400000.5, targetMjds, states, ephPlanet,
					                  planetRadiusKm(cfg.planet), *eph, nonGrav ? *nonGrav : propagator.nonGrav(),
					                  NBodyPropagator::SurfaceShape::Sphere, 0.0, control->cancel)
					            : propagator.propagateUntilSurface(r0Au, v0AuPerDay, epochTdb, targets, states,
					                                               ephPlanet, planetRadiusKm(cfg.planet), *eph,
					                                               nonGrav ? *nonGrav : propagator.nonGrav(),
					                                               NBodyPropagator::SurfaceShape::Sphere);
					if (control && control->cancel.stop_requested())
						throw std::runtime_error("Close-approach search cancelled");
					if (control && !terminal.hitSurface && states.size() != targets.size())
						throw std::runtime_error("Incomplete close-approach propagation");
					for (std::size_t i = 0; i < states.size(); ++i)
						appendSample(targets[i], states[i].R, states[i].V);
					if (terminal.hitSurface) {
						appendSample(terminal.jd, terminal.R, terminal.V);
						result.surfaceEntries.push_back(terminal);
						audit::write("grid-surface-stop", [&](auto& log) {
							log << " jd=" << terminal.jd << " completed=" << states.size();
						});
					}
				};
				const auto split = std::lower_bound(grid.begin(), grid.end(), epochTdb);
				std::vector<double> backward(grid.begin(), split);
				std::reverse(backward.begin(), backward.end());
				runDirection(backward);
				std::reverse(result.grid.begin(), result.grid.end());
				std::reverse(result.samples.begin(), result.samples.end());
				runDirection(std::vector<double>(split, grid.end()));
			}

			catch (const std::exception& error) {
				if (control)
					throw;
				audit::write("grid-failed", [&](auto& log) { log << " error=" << std::quoted(error.what()); });
				return std::nullopt;
			} catch (...) {
				if (control)
					throw;
				audit::write("grid-failed", [](auto& log) { log << " error=unknown"; });
				return std::nullopt;
			}

			if (result.grid.empty())
				return std::nullopt;
			audit::write("grid-complete", [&](auto& log) {
				log << " count=" << result.grid.size() << " startJd=" << result.grid.front()
				    << " endJd=" << result.grid.back();
			});
			return result;
		}

		[[nodiscard]] double radialVelocityAuPerDay(const Sample& s) noexcept
		{
			const Vec3D dr = s.rTarget - s.rPlanet;
			const double r = dr.magnitude();
			if (!(r > 0.0))
				return std::numeric_limits<double>::quiet_NaN();

			return dr.dot(s.vTarget - s.vPlanet) / r;
		}

		[[nodiscard]] std::optional<Sample> evaluateSample(const NBodyPropagator& propagator, IEphemerisSource* eph,
		                                                   EphemerisPlanet ephPlanet, const Sample& anchor,
		                                                   double anchorJd, double targetJd,
		                                                   const NonGravModel* nonGrav,
		                                                   const FitRunControl* control = nullptr)
		{
			const std::vector<double> grid{targetJd};
			std::vector<NBodyPropagator::StateVel> state;

			try {
				if (control && control->cancel.stop_requested())
					throw std::runtime_error("Close-approach search cancelled");
				state = control ? propagator.propagateAllWithVelocityMjd(
				                      anchor.rTarget, anchor.vTarget, anchorJd - 2400000.5, {targetJd - 2400000.5},
				                      nonGrav ? *nonGrav : propagator.nonGrav(), control->cancel)
				        : nonGrav ? propagator.propagateAllWithVelocity(anchor.rTarget, anchor.vTarget, anchorJd, grid,
				                                                        *nonGrav)
				                  : propagator.propagateAllWithVelocity(anchor.rTarget, anchor.vTarget, anchorJd, grid);
			} catch (...) {
				if (control)
					throw;
				return std::nullopt;
			}

			if (state.size() != 1)
				return std::nullopt;

			EphemerisRecord planetState;
			try {
				planetState = eph->queryJD(ephPlanet, targetJd);
			} catch (...) {
				if (control)
					throw;
				return std::nullopt;
			}

			Sample s;
			s.rTarget = state.front().R;
			s.vTarget = state.front().V;
			s.rPlanet = planetState.R;
			s.vPlanet = planetState.V;
			const Vec3D dr = s.rTarget - s.rPlanet;
			s.d2 = dr.dot(dr);
			return s;
		}

		void assignEventFromSample(CloseApproachEvent& ev, const Sample& s, double jd, Planet planet) noexcept
		{
			ev.jd = jd;
			ev.distAu = std::sqrt(s.d2);
			ev.planet = planet;
			ev.rTarget = s.rTarget;
			ev.vTarget = s.vTarget;
			ev.rPlanet = s.rPlanet;
			ev.vPlanet = s.vPlanet;

			const Vec3D vRel = ev.vTarget - ev.vPlanet;
			ev.vInfKmS = closeApproachVInfinityKmS(vRel, ev.distAu, planet);
		}

		[[nodiscard]] std::size_t encounterEntryIndex(const std::vector<Sample>& samples,
		                                              std::size_t closestApproachBracket, double thresholdAu2) noexcept
		{
			std::size_t entry = closestApproachBracket;
			while (entry > 0 && samples[entry - 1].d2 <= thresholdAu2)
				--entry;
			return entry;
		}

		[[nodiscard]] bool refineEventByRadialVelocity(const NBodyPropagator& propagator, IEphemerisSource* eph,
		                                               EphemerisPlanet ephPlanet, Planet planet,
		                                               const Sample& leftSample, double leftJd,
		                                               const Sample& rightSample, double rightJd,
		                                               const NonGravModel* nonGrav, CloseApproachEvent& event,
		                                               const FitRunControl* control = nullptr)
		{
			constexpr int kMaxIterations = 50;
			constexpr double kTimeToleranceDays = 1.0e-10;
			constexpr double kRdotToleranceAuPerDay = 1.0e-13;

			Sample s1 = leftSample;
			Sample s2 = rightSample;
			double t1 = leftJd;
			double t2 = rightJd;
			double rdot1 = radialVelocityAuPerDay(s1);
			double rdot2 = radialVelocityAuPerDay(s2);

			if (!std::isfinite(rdot1) || !std::isfinite(rdot2))
				return false;

			const double direction = (t2 > t1) ? 1.0 : -1.0;
			if (std::abs(rdot1) <= kRdotToleranceAuPerDay) {
				assignEventFromSample(event, s1, t1, planet);
				return true;
			}
			if (std::abs(rdot2) <= kRdotToleranceAuPerDay) {
				assignEventFromSample(event, s2, t2, planet);
				return true;
			}
			if (!(rdot1 * direction < 0.0 && rdot2 * direction > 0.0))
				return false;

			std::optional<Sample> lastSample;
			double lastJd = t2;

			for (int it = 0; it < kMaxIterations; ++it) {
				const double denom = rdot1 - rdot2;
				if (!control && std::abs(denom) <= std::numeric_limits<double>::epsilon())
					break;

				const double dt = -rdot2 * (t1 - t2) / denom;
				double t0 = t2 + dt;
				const double lo = std::min(t1, t2);
				const double hi = std::max(t1, t2);
				if (control)
					t0 = lo + (hi - lo) * 0.5;
				if (!(t0 >= lo && t0 <= hi))
					break;

				const std::optional<Sample> s0 =
				    evaluateSample(propagator, eph, ephPlanet, s2, t2, t0, nonGrav, control);
				if (!s0)
					break;

				const double rdot0 = radialVelocityAuPerDay(*s0);
				if (!std::isfinite(rdot0))
					break;

				lastSample = *s0;
				lastJd = t0;

				if (std::abs(rdot0) <= kRdotToleranceAuPerDay || std::abs(t2 - t1) <= kTimeToleranceDays ||
				    (control ? hi - lo <= 1e-9 : std::abs(dt) <= kTimeToleranceDays)) {
					assignEventFromSample(event, *s0, t0, planet);
					return true;
				}

				if (rdot0 * rdot1 < 0.0) {
					s2 = *s0;
					t2 = t0;
					rdot2 = rdot0;
				} else {
					s1 = *s0;
					t1 = t0;
					rdot1 = rdot0;
				}
			}

			if (control)
				throw std::runtime_error("Close-approach minimum refinement did not converge: " +
				                         std::string(longLabel(planet)) + " at JD " + std::to_string(leftJd));
			if (!lastSample)
				return false;

			assignEventFromSample(event, *lastSample, lastJd, planet);
			return true;
		}

		[[nodiscard]] std::vector<CloseApproachEvent> extractEvents(
		    const NBodyPropagator& propagator, IEphemerisSource* eph, const SampleGrid& sg,
		    const CloseApproachConfig& cfg, const NonGravModel* nonGrav,
		    std::vector<CloseApproachPassage>* passages = nullptr, const FitRunControl* control = nullptr)
		{
			std::vector<CloseApproachEvent> out;
			const std::vector<double>& grid = sg.grid;
			const std::vector<Sample>& samples = sg.samples;
			const double thresholdAu2 = cfg.thresholdAu * cfg.thresholdAu;
			const double encounterThresholdAu = cfg.encounterRadiusAu;
			const double encounterThresholdAu2 = encounterThresholdAu * encounterThresholdAu;
			const EphemerisPlanet ephPlanet = toEphemerisPlanet(cfg.planet);

			std::vector<int> encounterIndices(samples.size(), -1);
			int nextEncounterIndex = 0;
			bool insideEncounter = false;
			int currentEncounterIndex = -1;
			for (std::size_t i = 0; i < samples.size(); ++i) {
				if (samples[i].d2 <= encounterThresholdAu2) {
					if (!insideEncounter) {
						currentEncounterIndex = nextEncounterIndex++;
						insideEncounter = true;
						if (passages)
							passages->push_back({currentEncounterIndex, grid[i], grid[i]});
					}
					encounterIndices[i] = currentEncounterIndex;
					if (passages)
						passages->back().lastInsideJd = grid[i];
				} else {
					insideEncounter = false;
					currentEncounterIndex = -1;
				}
			}

			double prevRdot = radialVelocityAuPerDay(samples.front());
			for (std::size_t i = 1; i < samples.size(); ++i) {
				if (control && control->cancel.stop_requested())
					throw std::runtime_error("Close-approach search cancelled");
				const double currRdot = radialVelocityAuPerDay(samples[i]);
				if (!std::isfinite(prevRdot) || !std::isfinite(currRdot)) {
					prevRdot = currRdot;
					continue;
				}

				const bool inboundToOutbound =
				    control ? (prevRdot <= 0.0 && currRdot >= 0.0 && (prevRdot < 0.0 || currRdot > 0.0))
				            : (prevRdot < 0.0 && currRdot > 0.0);
				if (!inboundToOutbound) {
					prevRdot = currRdot;
					continue;
				}

				if (!control && std::min(samples[i - 1].d2, samples[i].d2) > thresholdAu2 * 4.0) {
					prevRdot = currRdot;
					continue;
				}

				CloseApproachEvent ev;
				if (refineEventByRadialVelocity(propagator, eph, ephPlanet, cfg.planet, samples[i - 1], grid[i - 1],
				                                samples[i], grid[i], nonGrav, ev, control) &&
				    ev.distAu < cfg.thresholdAu) {
					const std::size_t entry = encounterEntryIndex(samples, i - 1, encounterThresholdAu2);
					updateLocalMoid(ev, samples[entry], grid[entry]);
					if (encounterIndices[i - 1] >= 0) {
						ev.encounterIndex = encounterIndices[i - 1];
					} else if (encounterIndices[i] >= 0) {
						ev.encounterIndex = encounterIndices[i];
					} else {
						ev.encounterIndex = nextEncounterIndex++;
						if (passages)
							passages->push_back({ev.encounterIndex, ev.jd, ev.jd});
					}
					if (!control || out.empty() || std::abs(out.back().jd - ev.jd) > 1e-8)
						out.push_back(ev);
					audit::write("ca-record", [&](auto& log) {
						const Vec3D dr = ev.rTarget - ev.rPlanet;
						const Vec3D dv = ev.vTarget - ev.vPlanet;
						log << " encounter=" << ev.encounterIndex << " jd=" << ev.jd << " distanceAu=" << ev.distAu
						    << " rdot=" << dr.dot(dv) / dr.magnitude() << " rx=" << dr.x << " ry=" << dr.y
						    << " rz=" << dr.z << " vx=" << dv.x << " vy=" << dv.y << " vz=" << dv.z;
					});
				}

				prevRdot = currRdot;
			}

			for (const auto& surface : sg.surfaceEntries) {
				const auto p = eph->queryJD(ephPlanet, surface.jd);
				const auto mapping = mapInteriorEncounterToPericenter(surface.R - p.R, surface.V - p.V, cfg.planet);
				if (!mapping.valid)
					continue;
				const double jd = surface.jd + mapping.deltaDays;
				const auto periPlanet = eph->queryJD(ephPlanet, jd);
				const Sample peri{periPlanet.R + mapping.relativePosition, periPlanet.V + mapping.relativeVelocity,
				                  periPlanet.R, periPlanet.V, mapping.relativePosition.dot(mapping.relativePosition)};
				CloseApproachEvent ev;
				assignEventFromSample(ev, peri, jd, cfg.planet);
				if (ev.distAu >= cfg.thresholdAu)
					continue;
				const auto pos = std::lower_bound(grid.begin(), grid.end(), surface.jd);
				const std::size_t index = static_cast<std::size_t>(pos - grid.begin());
				const std::size_t entry = encounterEntryIndex(samples, index, encounterThresholdAu2);
				updateLocalMoid(ev, samples[entry], grid[entry]);
				if (encounterIndices[index] >= 0) {
					ev.encounterIndex = encounterIndices[index];
				} else {
					ev.encounterIndex = nextEncounterIndex++;
					if (passages)
						passages->push_back({ev.encounterIndex, surface.jd, surface.jd});
				}
				out.push_back(ev);
				audit::write("ca-record", [&](auto& log) {
					const Vec3D dr = ev.rTarget - ev.rPlanet;
					const Vec3D dv = ev.vTarget - ev.vPlanet;
					log << " encounter=" << ev.encounterIndex << " jd=" << ev.jd << " distanceAu=" << ev.distAu
					    << " rdot=" << dr.dot(dv) / dr.magnitude() << " rx=" << dr.x << " ry=" << dr.y << " rz=" << dr.z
					    << " vx=" << dv.x << " vy=" << dv.y << " vz=" << dv.z;
				});
			}
			std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.jd < b.jd; });
			return out;
		}
	}  // namespace

	std::vector<CloseApproachEvent> findCloseApproachesChecked(const Vec3D& r, const Vec3D& v, double epochJdTdb,
	                                                           const NBodyPropagator& propagator,
	                                                           const NonGravModel& ng,
	                                                           const CloseApproachConfig& config, IEphemerisSource* eph,
	                                                           const FitRunControl& control)
	{
		const auto grid = buildSampleGridFromState(r, v, epochJdTdb, propagator, &ng, config, eph, true, &control);
		if (!grid)
			throw std::runtime_error("Close-approach sampling failed");
		return extractEvents(propagator, eph, *grid, config, &ng, nullptr, &control);
	}

	std::vector<CloseApproachEvent> findCloseApproaches(const Orbit& orbit, const NBodyPropagator& propagator,
	                                                    const CloseApproachConfig& cfg, IEphemerisSource* eph)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		const std::optional<SampleGrid> sg = buildSampleGrid(orbit, propagator, cfg, eph);
		if (!sg)
			return {};

		return extractEvents(propagator, eph, *sg, cfg, nullptr);
	}

	std::vector<CloseApproachEvent> findCloseApproachesFromState(const Vec3D& r0Au, const Vec3D& v0AuPerDay,
	                                                             double epochTdb, const NBodyPropagator& propagator,
	                                                             const CloseApproachConfig& cfg, IEphemerisSource* eph)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		const std::optional<SampleGrid> sg =
		    buildSampleGridFromState(r0Au, v0AuPerDay, epochTdb, propagator, nullptr, cfg, eph);
		if (!sg)
			return {};

		return extractEvents(propagator, eph, *sg, cfg, nullptr);
	}

	std::vector<CloseApproachEvent> findCloseApproachesFromState(const Vec3D& r0Au, const Vec3D& v0AuPerDay,
	                                                             double epochTdb, const NBodyPropagator& propagator,
	                                                             const NonGravModel& nonGrav,
	                                                             const CloseApproachConfig& cfg, IEphemerisSource* eph)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		const std::optional<SampleGrid> sg =
		    buildSampleGridFromState(r0Au, v0AuPerDay, epochTdb, propagator, &nonGrav, cfg, eph);
		if (!sg)
			return {};

		return extractEvents(propagator, eph, *sg, cfg, &nonGrav);
	}

	std::vector<CloseApproachEvent> findCloseApproachesFromState(const Vec3D& r0Au, const Vec3D& v0AuPerDay,
	                                                             double epochTdb, const NBodyPropagator& propagator,
	                                                             const NonGravModel& nonGrav,
	                                                             const CloseApproachConfig& cfg, IEphemerisSource* eph,
	                                                             std::vector<CloseApproachPassage>& passages)
	{
		passages.clear();
		if (!eph)
			eph = getDefaultEphemerisSource();
		const auto sg = buildSampleGridFromState(r0Au, v0AuPerDay, epochTdb, propagator, &nonGrav, cfg, eph, true);
		if (!sg)
			return {};
		return extractEvents(propagator, eph, *sg, cfg, &nonGrav, &passages);
	}

}  // namespace orbital_core::hazard
