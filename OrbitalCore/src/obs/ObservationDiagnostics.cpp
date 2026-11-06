#include "orbital_core/ObservationDiagnostics.h"

#include "orbital_core/Constants.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"
#include "earth_orientation/EarthRotation.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>
#include <array>
#include <string_view>

namespace orbital_core {
	namespace {
		bool eligible(const ObsNode& o)
		{
			return !o.radar && !o.isSatellite && !o.rovingObserver && !o.obsCode.empty() && std::isfinite(o.JD) &&
			       std::isfinite(o.submittedRaRad()) && std::isfinite(o.submittedDecRad()) &&
			       (!o.hasSplitUtc || std::isfinite(o.utcSecondsOfDay));
		}

		double separation(const ObsNode& a, const ObsNode& b)
		{
			const double dra = std::remainder(a.submittedRaRad() - b.submittedRaRad(), Constants::TWOPI);
			const double ddec = a.submittedDecRad() - b.submittedDecRad();
			const double h = std::pow(std::sin(ddec / 2), 2) + std::cos(a.submittedDecRad()) *
			                                                       std::cos(b.submittedDecRad()) *
			                                                       std::pow(std::sin(dra / 2), 2);
			return 2 * std::asin(std::sqrt(std::clamp(h, 0.0, 1.0))) / Constants::ARCSEC_RAD;
		}
	}  // namespace

	std::vector<ObservationTimeGroup> findCoincidentObservations(std::span<const ObsNode> observations)
	{
		struct TimedObservation {
			std::size_t index;
			double mjdTt;
		};

		std::vector<TimedObservation> order;
		order.reserve(observations.size());
		for (std::size_t i = 0; i < observations.size(); ++i) {
			if (eligible(observations[i])) {
				const auto& o = observations[i];
				const double sourceMjd =
				    o.hasSplitUtc ? o.utcMjdDay + o.utcSecondsOfDay / 86400.0 : time::mjdFromJd(o.JD);
				const double tt = time::toMjdTt(sourceMjd, o.timeScale);
				if (std::isfinite(tt))
					order.push_back({i, tt});
			}
		}
		std::sort(order.begin(), order.end(), [&](const auto& a, const auto& b) {
			return std::tuple{a.mjdTt, observations[a.index].obsCode, observations[a.index].note1, a.index} <
			       std::tuple{b.mjdTt, observations[b.index].obsCode, observations[b.index].note1, b.index};
		});

		std::vector<ObservationTimeGroup> result;
		for (std::size_t first = 0; first < order.size();) {
			std::size_t end = first + 1;
			while (end < order.size() &&
			       observations[order[first].index].obsCode == observations[order[end].index].obsCode &&
			       (order[end].mjdTt - order[first].mjdTt) * 86400.0 < 3.0)
				++end;
			if (end - first > 1) {
				ObservationTimeGroup group;
				group.observations.reserve(end - first);
				const auto& reference = observations[order[first].index];
				std::set<std::tuple<double, double, double>> coordinates;
				bool samePosition = false;
				double minRa = 0.0, maxRa = 0.0, minDec = 0.0, maxDec = 0.0;
				for (std::size_t i = first; i < end; ++i) {
					const auto& observation = observations[order[i].index];
					group.observations.push_back(observation.id);
					group.maxSeparationArcsec = std::max(group.maxSeparationArcsec, separation(reference, observation));
					const auto normalizedRa = std::remainder(observation.submittedRaRad(), Constants::TWOPI);
					if (!coordinates.emplace(order[i].mjdTt, normalizedRa, observation.submittedDecRad()).second)
						samePosition = true;
					const double dra =
					    std::remainder(observation.submittedRaRad() - reference.submittedRaRad(), Constants::TWOPI);
					const double ddec = observation.submittedDecRad() - reference.submittedDecRad();
					minRa = std::min(minRa, dra);
					maxRa = std::max(maxRa, dra);
					minDec = std::min(minDec, ddec);
					maxDec = std::max(maxDec, ddec);
				}
				if (samePosition)
					group.issues.push_back(ObservationIssue::CoincidentPosition);
				if (maxRa - minRa > 2.0 * Constants::ARCSEC_RAD || maxDec - minDec > 2.0 * Constants::ARCSEC_RAD)
					group.issues.push_back(ObservationIssue::NearbyObservationConflict);
				if (!group.issues.empty())
					result.push_back(std::move(group));
			}
			first = end;
		}
		return result;
	}

	std::vector<ObservationVisibilityDiagnostic> diagnoseObservationVisibility(std::span<const ObsNode> observations)
	{
		// Radio receivers that can observe in daylight. Match complete station codes.
		constexpr std::array<std::string_view, 4> daytimeStations{"ScT", "ATA", "256", "EKa"};
		std::vector<ObservationVisibilityDiagnostic> result;
		for (const auto& o : observations) {
			if (o.radar || o.isSatellite || o.obsCode.empty())
				continue;
			ObservationVisibilityDiagnostic item;
			item.observation = o.id;
			const double rho = std::hypot(o.cosPhi, o.sinPhi);
			if (std::isfinite(rho) && rho > 0.0 && std::isfinite(o.LST) && std::isfinite(o.submittedRaRad()) &&
			    std::isfinite(o.submittedDecRad())) {
				try {
					const double jdTt = time::jdFromMjd(time::fromMjdTdb(o.epochMjdTdb(), time::Scale::TT));
					const double e2 = Constants::EARTH_FLATTENING * (2.0 - Constants::EARTH_FLATTENING);
					const double latitude = std::atan2(o.sinPhi, o.cosPhi * (1.0 - e2));
					const Vec3D zenith(std::cos(latitude) * std::cos(o.LST), std::cos(latitude) * std::sin(o.LST),
					                   std::sin(latitude));
					const double ra = o.submittedRaRad(), dec = o.submittedDecRad();
					Vec3D target(std::cos(dec) * std::cos(ra), std::cos(dec) * std::sin(ra), std::sin(dec));
					Vec3D sun = Helpers::eclToEcv(o.obsVector);
					EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, target, 1);
					EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, sun, 1);
					const auto altitude = [&](const Vec3D& direction) -> std::optional<double> {
						const double length = direction.magnitude();
						if (!std::isfinite(length) || length <= 0.0)
							return {};
						return std::asin(std::clamp(direction.dot(zenith) / length, -1.0, 1.0)) * Constants::RAD_TO_DEG;
					};
					item.objectAltitudeDeg = altitude(target);
					item.sunAltitudeDeg = altitude(sun);
					if (item.objectAltitudeDeg && *item.objectAltitudeDeg < 0.0)
						item.issues.push_back(ObservationIssue::BelowHorizon);
					if (item.sunAltitudeDeg && *item.sunAltitudeDeg > 0.0 &&
					    std::ranges::find(daytimeStations, o.obsCode) == daytimeStations.end())
						item.issues.push_back(ObservationIssue::Daylight);
				} catch (const std::exception&) {
					// Missing Earth-orientation data makes altitude unavailable.
				}
			}
			if (!item.issues.empty())
				result.push_back(std::move(item));
		}
		return result;
	}
}  // namespace orbital_core
