#include "orbital_core/hazard/ObsMonteCarlo.h"
#include <random>
#include <map>

#include "orbital_core/Constants.h"
#include "orbital_core/FitSession.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/BPlane.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "obs/AstrometricSigma.h"
#include "hazard/GaussianSampler.h"
#include "NBodyPropagator.h"
#include "logger/Logger.h"

#include <Eigen/Cholesky>

namespace orbital_core::hazard {
	namespace {
		using Vec6 = Eigen::Matrix<double, 6, 1>;
		using Mat6 = Eigen::Matrix<double, 6, 6>;
		using VecX = Eigen::VectorXd;
		using MatX = Eigen::MatrixXd;

		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
		constexpr double kGaussK = Constants::GAUSS_K;
		constexpr double kEncounterHalfWindowDays = 1.0;

		void validateAltitude(double altitudeKm)
		{
			if (!std::isfinite(altitudeKm) || altitudeKm < 0.0)
				throw std::invalid_argument("MC altitude must be finite and nonnegative");
			if (altitudeKm > 0.0)
				LOG_INFO << "[ObsMC] boundary altitude=" << altitudeKm
				         << " km above reference ellipsoid; hits and P mean boundary entry, not ground impact";
		}

		[[nodiscard]] Vec6 nominalPhysicalState(const Orbit& orbit)
		{
			const Vec3D r = orbit.getR();
			const Vec3D v = orbit.getV() * kGaussK;
			Vec6 x;
			x << r.x, r.y, r.z, v.x, v.y, v.z;
			return x;
		}

		[[nodiscard]] VecX nominalPhysicalParameters(const Orbit& orbit, const NonGravModel& ng)
		{
			VecX x(6 + ng.numFitParams());
			x.head<6>() = nominalPhysicalState(orbit);

			for (int k = 0; k < ng.numFitParams(); ++k)
				x(6 + k) = ng.fitParameterValue(k);

			return x;
		}

		[[nodiscard]] NonGravModel nonGravFromPhysicalParameters(const NonGravModel& nominal, const VecX& params)
		{
			NonGravModel ng = nominal;
			const int n =
			    std::min<int>(ng.numFitParams(), static_cast<int>(std::max<Eigen::Index>(0, params.size() - 6)));

			for (int k = 0; k < n; ++k)
				ng.setFitParameterValue(k, params(6 + k));

			return ng;
		}

		[[nodiscard]] Mat6 stabilizedCovariance(const StateCov6x6& cov)
		{
			Mat6 out = cov;
			out = 0.5 * (out + out.transpose()).eval();

			for (int i = 0; i < 6; ++i) {
				const double diag = out(i, i);
				if (!std::isfinite(diag) || diag <= 0.0)
					out(i, i) = std::max(std::abs(diag), 1.0e-30);
			}

			return out;
		}

		[[nodiscard]] MatX stabilizedCovariance(const FitParameterCovariance& cov)
		{
			MatX out = cov;
			out = 0.5 * (out + out.transpose()).eval();

			for (Eigen::Index i = 0; i < out.rows(); ++i) {
				const double diag = out(i, i);
				if (!std::isfinite(diag) || diag <= 0.0)
					out(i, i) = std::max(std::abs(diag), 1.0e-30);
			}

			return out;
		}

		[[nodiscard]] Mat6 covarianceFactor(const StateCov6x6& cov)
		{
			Mat6 stable = stabilizedCovariance(cov);
			for (int attempt = 0; attempt < 8; ++attempt) {
				Eigen::LLT<Mat6> llt(stable);
				if (llt.info() == Eigen::Success)
					return llt.matrixL();

				const double jitter = std::pow(10.0, -30.0 + attempt * 2.0);
				stable.diagonal().array() += jitter;
			}

			return Mat6::Zero();
		}

		[[nodiscard]] MatX covarianceFactor(const FitParameterCovariance& cov)
		{
			MatX stable = stabilizedCovariance(cov);
			for (int attempt = 0; attempt < 8; ++attempt) {
				Eigen::LLT<MatX> llt(stable);
				if (llt.info() == Eigen::Success)
					return llt.matrixL();

				const double jitter = std::pow(10.0, -30.0 + attempt * 2.0);
				stable.diagonal().array() += jitter;
			}

			return MatX::Zero(cov.rows(), cov.cols());
		}

		[[nodiscard]] bool usableFactor(const Mat6& factor) noexcept
		{
			return factor.array().isFinite().all() && !factor.isZero(0.0);
		}

		[[nodiscard]] bool usableFactor(const MatX& factor) noexcept
		{
			return factor.array().isFinite().all() && !factor.isZero(0.0);
		}

		void recomputeDirectionVector(ObsNode& o) noexcept
		{
			const Vec3D eqv(std::cos(o.DEC) * std::cos(o.RA), std::cos(o.DEC) * std::sin(o.RA), std::sin(o.DEC));

			o.vector = Convertors::ecvToEcl(eqv, Constants::EARTH_INCL);
		}
	}  // namespace

	namespace {
		void logCorridorStatistics(const ObsMcResult& out, std::span<const CloseApproachEvent> events)
		{
			if (out.corridor.empty())
				return;

			const auto pct = [](std::vector<double> v, double p) {
				if (v.empty())
					return std::numeric_limits<double>::quiet_NaN();
				std::sort(v.begin(), v.end());
				const double x = p * static_cast<double>(v.size() - 1);
				const std::size_t i = static_cast<std::size_t>(x);
				const double f = x - static_cast<double>(i);
				return i + 1 < v.size() ? v[i] * (1.0 - f) + v[i + 1] * f : v.back();
			};
			const auto line = [&](const char* name, std::vector<double> v, const char* unit) {
				if (v.empty())
					return;
				LOG_INFO << "[ObsMC-CORRIDOR]     " << name << " [" << unit << "]"
				         << "  min=" << pct(v, 0.0) << "  p05=" << pct(v, 0.05) << "  p25=" << pct(v, 0.25)
				         << "  med=" << pct(v, 0.50) << "  p75=" << pct(v, 0.75) << "  p95=" << pct(v, 0.95)
				         << "  max=" << pct(v, 1.0);
			};

			for (std::size_t e = 0; e < events.size(); ++e) {
				std::vector<double> tMin, lon, lat, q, vEnt, vInf;
				const double tNominal = events[e].jd;
				const double rPlanetKm = planetRadiusKm(events[e].planet);
				const double muKm3 = planetMuKm3PerS2(events[e].planet);
				const double vEsc = (rPlanetKm > 0.0 && muKm3 > 0.0) ? std::sqrt(2.0 * muKm3 / rPlanetKm) : 0.0;
				for (const auto& c : out.corridor) {
					if (static_cast<std::size_t>(c.eventIndex) != e)
						continue;
					tMin.push_back((c.tcaJd - tNominal) * 1440.0);
					lon.push_back(c.longitudeDeg);
					lat.push_back(c.latitudeDeg);
					q.push_back(c.missKm);
					vEnt.push_back(c.velocityKmS);
					const double d = c.velocityKmS * c.velocityKmS - vEsc * vEsc;
					vInf.push_back(d > 0.0 ? std::sqrt(d) : 0.0);
				}
				if (tMin.empty())
					continue;
				const double tSpanMin = pct(tMin, 1.0) - pct(tMin, 0.0);
				const double lonSpan = pct(lon, 1.0) - pct(lon, 0.0);
				LOG_INFO << "[ObsMC-CORRIDOR] event[" << e << "]"
				         << " points=" << tMin.size() << "  tNominalJd=" << std::setprecision(12) << tNominal
				         << std::setprecision(6) << "  vEsc=" << vEsc << " km/s"
				         << "  rPlanet=" << rPlanetKm << " km";
				LOG_INFO << "[ObsMC-CORRIDOR]   span: time=" << tSpanMin << " min"
				         << "  longitude=" << lonSpan << " deg";
				line("t-tNominal", tMin, "min");
				line("longitude ", lon, "deg");
				line("latitude  ", lat, "deg");
				line("pericentre", q, "km");
				line("entrySpeed", vEnt, "km/s");
				line("vInf      ", vInf, "km/s");
			}
		}
	}  // namespace

	namespace {
		struct CloneTrajectory {
			Vec3D r;
			Vec3D v;
			double epochJd;
			const NonGravModel& nonGrav;
		};

		struct EncounterSample {
			double missKm{std::numeric_limits<double>::quiet_NaN()};
			double conicPericenterKm{std::numeric_limits<double>::quiet_NaN()};
			std::optional<ImpactCloneInfo> impact;
		};

		EncounterSample describeEncounter(const CloseApproachEvent& event, const NBodyPropagator::StateVel* atEvent,
		                                  const NBodyPropagator::PropagationWithSurfaceStop& terminal,
		                                  IEphemerisSource& eph, double altitudeKm)
		{
			const auto body = toEphemerisPlanet(event.planet);
			EncounterSample result;
			if (terminal.hitSurface && terminal.jd < event.jd - kEncounterHalfWindowDays)
				return result;

			const double sampleJd = atEvent ? event.jd : terminal.jd;
			const auto sampleR = atEvent ? atEvent->R : terminal.R;
			const auto sampleV = atEvent ? atEvent->V : terminal.V;
			const auto planet = eph.queryJD(body, sampleJd);
			const Vec3D r = (sampleR - planet.R) * kAuKm;
			const Vec3D v = (sampleV - planet.V) * (kAuKm / 86400.0);
			const double rMag = r.magnitude();
			const double mu = planetMuKm3PerS2(event.planet);
			if (!std::isfinite(rMag) || !std::isfinite(v.dot(v)) || !(rMag > 0.0))
				throw std::runtime_error("Invalid MC encounter state");
			const Vec3D h = r.cross(v);
			const double h2 = h.dot(h);
			if (mu > 0.0) {
				const Vec3D e = v.cross(h) * (1.0 / mu) - r * (1.0 / rMag);
				result.conicPericenterKm = h2 / (mu * (1.0 + e.magnitude()));
			}

			const auto basis = computeBPlaneBasis(event);
			const double gamma = gravitationalFocusing(event.vInfKmS, event.planet);
			result.missKm = rMag / std::max(gamma, 1.0);
			if (basis.valid) {
				const double xi = r.dot(basis.xi) / gamma;
				const double zeta = r.dot(basis.zeta) / gamma;
				result.missKm = std::sqrt(xi * xi + zeta * zeta);
			}

			if (terminal.hitSurface && terminal.jd <= event.jd + kEncounterHalfWindowDays) {
				const auto atCrossing =
				    altitudeKm > 0.0 ? eph.queryMjd(body, terminal.mjd) : eph.queryJD(body, terminal.jd);
				const Vec3D relativeR = (terminal.R - atCrossing.R) * kAuKm;
				const Vec3D relativeV = terminal.V - atCrossing.V;
				if (!std::isfinite(relativeR.dot(relativeR)) || !std::isfinite(relativeV.dot(relativeV)))
					throw std::runtime_error("Invalid MC surface state");
				ImpactCloneInfo info;
				info.tcaJd = terminal.jd;
				info.planet = event.planet;
				info.missKm = result.conicPericenterKm;
				info.velocityKmS = relativeV.magnitude() * (kAuKm / 86400.0);
				info.rBodyFixedKm = icrfToBodyFixedMjd(event.planet, relativeR, time::mjdFromJd(terminal.jd), &eph);
				info.epochMjdTdb = terminal.mjd;
				info.epochMjdUtc = time::fromMjdTdb(terminal.mjd, time::Scale::UTC);
				const double beforeJd = terminal.jd - 1.0 / 86400.0;
				const double afterJd = terminal.jd + 1.0 / 86400.0;
				const Vec3D frameRate = (icrfToBodyFixedMjd(event.planet, relativeR, time::mjdFromJd(afterJd), &eph) -
				                         icrfToBodyFixedMjd(event.planet, relativeR, time::mjdFromJd(beforeJd), &eph)) *
				                        (1.0 / ((afterJd - beforeJd) * 86400.0));
				info.vBodyFixedKmS = icrfToBodyFixedMjd(event.planet, relativeV * (kAuKm / 86400.0),
				                                        time::mjdFromJd(terminal.jd), &eph) +
				                     frameRate;
				const auto ll = latLonFromBodyFixed(event.planet, info.rBodyFixedKm);
				info.latitudeDeg = ll.latitudeDeg;
				info.longitudeDeg = ll.longitudeDeg;
				info.bodyFixedValid = ll.valid;
				result.impact = info;
			}
			return result;
		}

		std::vector<EncounterSample> sampleEncounters(const CloneTrajectory& clone,
		                                              std::span<const CloseApproachEvent> events,
		                                              const NBodyPropagator& propagator, IEphemerisSource* eph,
		                                              double altitudeKm)
		{
			if (!eph || !eph->isAvailable())
				throw std::runtime_error("MC requires an available ephemeris source");
			std::map<Planet, std::vector<std::size_t>> byBody;
			for (std::size_t i = 0; i < events.size(); ++i) {
				if (!std::isfinite(events[i].jd) || events[i].jd < clone.epochJd)
					throw std::runtime_error("MC requires finite forward encounter epochs");
				byBody[events[i].planet].push_back(i);
			}
			std::vector<EncounterSample> result(events.size());
			for (const auto& [planet, indices] : byBody) {
				std::vector<double> targets;
				targets.reserve(indices.size() + 1);
				for (const auto index : indices)
					targets.push_back(events[index].jd);
				std::sort(targets.begin(), targets.end());
				targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
				targets.push_back(targets.back() + kEncounterHalfWindowDays);
				std::vector<NBodyPropagator::StateVel> samples;
				const auto terminal = propagator.propagateUntilSurface(
				    clone.r, clone.v, clone.epochJd, targets, samples, toEphemerisPlanet(planet),
				    planetRadiusKm(planet), *eph, clone.nonGrav,
				    altitudeKm > 0.0 ? NBodyPropagator::SurfaceShape::ReferenceEllipsoid
				                     : NBodyPropagator::SurfaceShape::Sphere,
				    altitudeKm);
				if (altitudeKm > 0.0 && terminal.hitSurface && terminal.jd == clone.epochJd)
					throw std::runtime_error("MC clone starts inside the requested altitude boundary");
				if (!terminal.hitSurface && (terminal.jd != targets.back() || samples.size() != targets.size()))
					throw std::runtime_error("Incomplete MC surface propagation");
				for (const auto index : indices) {
					const auto sampleIndex = static_cast<std::size_t>(
					    std::lower_bound(targets.begin(), targets.end(), events[index].jd) - targets.begin());
					const auto* atEvent = sampleIndex < samples.size() ? &samples[sampleIndex] : nullptr;
					result[index] = describeEncounter(events[index], atEvent, terminal, *eph, altitudeKm);
				}
			}
			return result;
		}

		void assignSurfacePoint(ObsMcResult& out, const Vec3D& R0, const Vec3D& V0, double epoch,
		                        const NonGravModel* ng, const NBodyPropagator& propagator,
		                        std::span<const CloseApproachEvent> events, IEphemerisSource* eph, double altitudeKm)
		{
			if (!eph || !eph->isAvailable())
				return;
			const CloneTrajectory clone{R0, V0, epoch, ng ? *ng : propagator.nonGrav()};
			try {
				for (const auto& sample : sampleEncounters(clone, events, propagator, eph, altitudeKm)) {
					if (!sample.impact || !sample.impact->bodyFixedValid)
						continue;
					out.centroidLatDeg = sample.impact->latitudeDeg;
					out.centroidLonDeg = sample.impact->longitudeDeg;
					out.centroidTcaJd = sample.impact->tcaJd;
					return;
				}
			} catch (...) {
			}
		}

		void assignNominalSurfacePoint(ObsMcResult& out, const Orbit& nominalOrbit, const NBodyPropagator& propagator,
		                               std::span<const CloseApproachEvent> events, IEphemerisSource* eph,
		                               double altitudeKm)
		{
			const Vec3D R0 = nominalOrbit.getR();
			const Vec3D V0 = nominalOrbit.getV() * kGaussK;
			assignSurfacePoint(out, R0, V0, nominalOrbit.getEpoch(), nullptr, propagator, events, eph, altitudeKm);
		}
	}  // namespace

	ObsMcResult runStateCovarianceMonteCarlo(const Orbit& nominalOrbit, const StateCov6x6& stateCovariance,
	                                         const NBodyPropagator& propagator,
	                                         std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg,
	                                         IEphemerisSource* eph, std::function<bool(int, int)> progressCb)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		validateAltitude(cfg.altitudeKm);
		ObsMcResult out;
		out.earthFrameModel = earthBodyFixedFrameModel();
		out.altitudeKm = cfg.altitudeKm;
		const int N = std::max(0, cfg.numClones);
		const std::size_t Nev = events.size();

		out.numAttempted = N;
		out.hitsPerEvent.assign(Nev, 0);
		out.minMissPerEventKm.assign(Nev, std::numeric_limits<double>::max());
		out.medianMissPerEventKm.assign(Nev, 0.0);

		if (N == 0 || events.empty() || stateCovariance.isZero(0.0))
			return out;

		const Vec6 xNom = nominalPhysicalState(nominalOrbit);
		Mat6 factor = covarianceFactor(stateCovariance);
		if (!usableFactor(factor))
			return out;
		factor *= (cfg.noiseSigmaMult > 0.0) ? cfg.noiseSigmaMult : 1.0;

		{
			const Vec3D R = nominalOrbit.getR();
			const Vec3D V = nominalOrbit.getV();
			LOG_INFO << "[ObsMC] start  rngSeed=" << cfg.rngSeed << "  N=" << N << "  sigmaMult=" << cfg.noiseSigmaMult
			         << "  Nev=" << events.size() << "  epoch=" << std::fixed << std::setprecision(9)
			         << nominalOrbit.getEpoch() << "  R=(" << R.x << "," << R.y << "," << R.z << ")"
			         << "  V=(" << V.x << "," << V.y << "," << V.z << ")"
			         << "  technique=cov-MC";

			std::ostringstream diag;
			diag << std::scientific << std::setprecision(3);
			for (int i = 0; i < 6; ++i)
				diag << " " << std::sqrt(std::max(0.0, stateCovariance(i, i)));
			LOG_INFO << "[ObsMC] cov 1sigma diag (R AU, V AU/d):" << diag.str();
		}

		std::vector<double> eventGamma(Nev);
		for (std::size_t e = 0; e < Nev; ++e) {
			eventGamma[e] = gravitationalFocusing(events[e].vInfKmS, events[e].planet);

			const double rPl = planetRadiusKm(events[e].planet);
			LOG_INFO << "[ObsMC] event[" << e << "] planet=" << static_cast<int>(events[e].planet)
			         << "  v_inf=" << std::fixed << std::setprecision(2) << events[e].vInfKmS << " km/s"
			         << "  gamma=" << std::setprecision(4) << eventGamma[e] << "  R=" << std::setprecision(0) << rPl
			         << " km"
			         << "  surface-radius=" << rPl << " km"
			         << "  (gammaR=" << eventGamma[e] * rPl << " km, info only)";
		}

		std::vector<std::vector<double>> missAccum(Nev);
		std::vector<std::array<double, 6>> cloneStates;

		std::atomic<int> numPropagated{0};
		std::atomic<int> numFailed{0};
		std::atomic<int> numAnyImpact{0};
		std::atomic<int> doneCount{0};
		std::atomic<bool> cancelFlag{false};

		std::atomic<int> totalImpacts{0};
		std::atomic<int> totalCorridorEntries{0};

		const std::size_t slotCount = static_cast<std::size_t>(N) * Nev;
		std::vector<double> missFlat(slotCount, std::numeric_limits<double>::quiet_NaN());
		std::vector<char> impactFlat(slotCount, 0);
		std::vector<char> cloneConverged(static_cast<std::size_t>(N), 0);
		std::vector<std::vector<ImpactCloneInfo>> corridorByClone(static_cast<std::size_t>(N));

#pragma omp parallel for schedule(dynamic, 16)
		for (int k = 0; k < N; ++k) {
			if (cancelFlag.load(std::memory_order_relaxed))
				continue;

			Vec6 z;
			sampling::fillStandardNormal(z, cfg.rngSeed, k);
			const Vec6 x = xNom + factor * z;

			const Vec3D R0(x(0), x(1), x(2));
			const Vec3D V0(x(3), x(4), x(5));
			const double epoch = nominalOrbit.getEpoch();

			std::vector<EncounterSample> samples;
			try {
				samples =
				    sampleEncounters({R0, V0, epoch, propagator.nonGrav()}, events, propagator, eph, cfg.altitudeKm);
			} catch (...) {
				numFailed.fetch_add(1, std::memory_order_relaxed);
				const int done = doneCount.fetch_add(1, std::memory_order_relaxed) + 1;
				if (progressCb) {
#pragma omp critical(covmc_progress)
					{
						if (!progressCb(done, N))
							cancelFlag.store(true);
					}
				}
				continue;
			}

			numPropagated.fetch_add(1, std::memory_order_relaxed);
			cloneConverged[static_cast<std::size_t>(k)] = 1;

#pragma omp critical(covmc_state)
			{
				cloneStates.push_back({R0.x, R0.y, R0.z, V0.x, V0.y, V0.z});
			}

			bool anyHit = false;
			for (std::size_t e = 0; e < Nev; ++e) {
				auto& sample = samples[e];
				const std::size_t slot = static_cast<std::size_t>(k) * Nev + e;
				missFlat[slot] = sample.missKm;
				if (sample.impact) {
					totalImpacts.fetch_add(1, std::memory_order_relaxed);
					anyHit = true;
					impactFlat[slot] = 1;
					sample.impact->cloneIndex = k;
					sample.impact->eventIndex = static_cast<int>(e);
					corridorByClone[static_cast<std::size_t>(k)].push_back(std::move(*sample.impact));
				}
			}

			if (anyHit)
				numAnyImpact.fetch_add(1, std::memory_order_relaxed);

			const int done = doneCount.fetch_add(1, std::memory_order_relaxed) + 1;
			if (progressCb) {
#pragma omp critical(covmc_progress)
				{
					if (!progressCb(done, N))
						cancelFlag.store(true);
				}
			}
		}

		for (int k = 0; k < N; ++k) {
			if (!cloneConverged[static_cast<std::size_t>(k)])
				continue;
			for (std::size_t e = 0; e < Nev; ++e) {
				const std::size_t slot = static_cast<std::size_t>(k) * Nev + e;
				const double miss = missFlat[slot];
				if (std::isnan(miss))
					continue;
				missAccum[e].push_back(miss);
				if (miss < out.minMissPerEventKm[e])
					out.minMissPerEventKm[e] = miss;
				if (impactFlat[slot])
					++out.hitsPerEvent[e];
			}
			for (ImpactCloneInfo& info : corridorByClone[static_cast<std::size_t>(k)]) {
				out.corridor.push_back(std::move(info));
				totalCorridorEntries.fetch_add(1, std::memory_order_relaxed);
			}
		}

		out.numConverged = numPropagated.load();
		out.numFailedFit = numFailed.load();
		out.numAnyImpactClones = numAnyImpact.load();
		out.cancelled = cancelFlag.load();

		if (out.numConverged > 0) {
			const double n = static_cast<double>(out.numConverged);
			const double pAny = static_cast<double>(out.numAnyImpactClones) / n;
			out.aggregateImpactProb = pAny;
			const double binom = std::max(pAny * (1.0 - pAny), 3.0 / n);
			out.aggregateImpactSigma = std::sqrt(binom / n);
		}

		for (std::size_t e = 0; e < Nev; ++e) {
			auto& v = missAccum[e];
			if (v.empty())
				continue;

			std::sort(v.begin(), v.end());
			out.medianMissPerEventKm[e] = v[v.size() / 2];
		}

		for (auto& m : out.minMissPerEventKm) {
			if (m == std::numeric_limits<double>::max())
				m = std::numeric_limits<double>::quiet_NaN();
		}

		assignNominalSurfacePoint(out, nominalOrbit, propagator, events, eph, cfg.altitudeKm);

		const double convPct = (N > 0) ? 100.0 * static_cast<double>(out.numConverged) / static_cast<double>(N) : 0.0;

		LOG_INFO << "[ObsMC] N=" << N << "  converged=" << out.numConverged << " (" << std::fixed
		         << std::setprecision(0) << convPct << "%)"
		         << "  any-impact=" << out.numAnyImpactClones << "  P(any)=" << std::scientific << std::setprecision(3)
		         << out.aggregateImpactProb << " +/- " << out.aggregateImpactSigma;

		for (std::size_t e = 0; e < Nev; ++e) {
			LOG_INFO << "[ObsMC]   event[" << e << "] minMiss=" << std::fixed << std::setprecision(0)
			         << out.minMissPerEventKm[e] << " km"
			         << "  medianMiss=" << out.medianMissPerEventKm[e] << " km"
			         << "  hits=" << out.hitsPerEvent[e] << "/" << out.numConverged;
		}

		LOG_INFO << "[ObsMC] DIAG-XING summary: total impacts=" << totalImpacts.load()
		         << "  corridor entries added=" << totalCorridorEntries.load();

		logCorridorStatistics(out, events);

		if (cloneStates.size() >= 2) {
			std::array<double, 6> mean{};
			for (const auto& s : cloneStates) {
				for (int i = 0; i < 6; ++i)
					mean[i] += s[i];
			}
			for (int i = 0; i < 6; ++i)
				mean[i] /= static_cast<double>(cloneStates.size());

			std::array<double, 6> var{};
			for (const auto& s : cloneStates) {
				for (int i = 0; i < 6; ++i) {
					const double d = s[i] - mean[i];
					var[i] += d * d;
				}
			}

			const double denom = static_cast<double>(cloneStates.size() - 1);
			std::ostringstream sigStr;
			sigStr << std::scientific << std::setprecision(3);
			for (int i = 0; i < 6; ++i)
				sigStr << " " << std::sqrt(var[i] / denom);
			LOG_INFO << "[ObsMC]   empirical clone sigma (R AU, V AU/d):" << sigStr.str();
		}

		return out;
	}

	ObsMcResult runParameterCovarianceMonteCarloAtCenter(
	    const Orbit& nominalOrbit, const Eigen::VectorXd& parameterCenter,
	    const FitParameterCovariance& parameterCovariance, const NonGravModel& nominalNonGrav,
	    const NBodyPropagator& propagator, std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg,
	    IEphemerisSource* eph, std::function<bool(int, int)> progressCb)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		validateAltitude(cfg.altitudeKm);
		ObsMcResult out;
		out.earthFrameModel = earthBodyFixedFrameModel();
		out.altitudeKm = cfg.altitudeKm;
		const int N = std::max(0, cfg.numClones);
		const std::size_t Nev = events.size();
		const Eigen::Index dim = parameterCenter.size();
		const int expectedParams = 6 + nominalNonGrav.numFitParams();

		out.numAttempted = N;
		out.hitsPerEvent.assign(Nev, 0);
		out.minMissPerEventKm.assign(Nev, std::numeric_limits<double>::max());
		out.medianMissPerEventKm.assign(Nev, 0.0);

		if (N == 0 || events.empty() || dim < 6 || dim != expectedParams || parameterCovariance.rows() != dim ||
		    parameterCovariance.cols() != dim || parameterCovariance.isZero(0.0)) {
			return out;
		}

		MatX factor = covarianceFactor(parameterCovariance);
		if (!usableFactor(factor))
			return out;
		factor *= (cfg.noiseSigmaMult > 0.0) ? cfg.noiseSigmaMult : 1.0;

		{
			const Vec3D R(parameterCenter(0), parameterCenter(1), parameterCenter(2));
			const Vec3D V(parameterCenter(3) / kGaussK, parameterCenter(4) / kGaussK, parameterCenter(5) / kGaussK);
			LOG_INFO << "[ObsMC] start  rngSeed=" << cfg.rngSeed << "  N=" << N << "  sigmaMult=" << cfg.noiseSigmaMult
			         << "  Nev=" << events.size() << "  epoch=" << std::fixed << std::setprecision(9)
			         << nominalOrbit.getEpoch() << "  R=(" << R.x << "," << R.y << "," << R.z << ")"
			         << "  V=(" << V.x << "," << V.y << "," << V.z << ")"
			         << "  params=" << dim << "  technique=param-cov-MC";

			std::ostringstream diag;
			diag << std::scientific << std::setprecision(3);
			for (int i = 0; i < 6; ++i)
				diag << " " << std::sqrt(std::max(0.0, parameterCovariance(i, i)));
			LOG_INFO << "[ObsMC] cov 1sigma diag (R AU, V AU/d):" << diag.str();

			if (dim > 6) {
				std::ostringstream ngDiag;
				ngDiag << std::scientific << std::setprecision(3);
				for (Eigen::Index i = 6; i < dim; ++i)
					ngDiag << " " << std::sqrt(std::max(0.0, parameterCovariance(i, i)));
				LOG_INFO << "[ObsMC] cov 1sigma diag (non-grav params):" << ngDiag.str();
			}
		}

		std::vector<double> eventGamma(Nev);
		for (std::size_t e = 0; e < Nev; ++e) {
			eventGamma[e] = gravitationalFocusing(events[e].vInfKmS, events[e].planet);

			const double rPl = planetRadiusKm(events[e].planet);
			LOG_INFO << "[ObsMC] event[" << e << "] planet=" << static_cast<int>(events[e].planet)
			         << "  v_inf=" << std::fixed << std::setprecision(2) << events[e].vInfKmS << " km/s"
			         << "  gamma=" << std::setprecision(4) << eventGamma[e] << "  R=" << std::setprecision(0) << rPl
			         << " km"
			         << "  surface-radius=" << rPl << " km"
			         << "  (gammaR=" << eventGamma[e] * rPl << " km, info only)";
		}

		std::vector<std::vector<double>> missAccum(Nev);
		std::vector<std::array<double, 6>> cloneStates;

		std::atomic<int> numPropagated{0};
		std::atomic<int> numFailed{0};
		std::atomic<int> numAnyImpact{0};
		std::atomic<int> doneCount{0};
		std::atomic<bool> cancelFlag{false};

		std::atomic<int> totalImpacts{0};
		std::atomic<int> totalCorridorEntries{0};

		const std::size_t slotCount = static_cast<std::size_t>(N) * Nev;
		std::vector<double> missFlat(slotCount, std::numeric_limits<double>::quiet_NaN());
		std::vector<char> impactFlat(slotCount, 0);
		std::vector<char> cloneConverged(static_cast<std::size_t>(N), 0);
		std::vector<std::vector<ImpactCloneInfo>> corridorByClone(static_cast<std::size_t>(N));

#pragma omp parallel for schedule(dynamic, 16)
		for (int k = 0; k < N; ++k) {
			if (cancelFlag.load(std::memory_order_relaxed))
				continue;

			VecX z(dim);
			sampling::fillStandardNormal(z, cfg.rngSeed, k);
			const VecX x = parameterCenter + factor * z;

			const Vec3D R0(x(0), x(1), x(2));
			const Vec3D V0(x(3), x(4), x(5));
			const double epoch = nominalOrbit.getEpoch();
			const NonGravModel ng = nonGravFromPhysicalParameters(nominalNonGrav, x);

			std::vector<EncounterSample> samples;
			try {
				samples = sampleEncounters({R0, V0, epoch, ng}, events, propagator, eph, cfg.altitudeKm);
			} catch (...) {
				numFailed.fetch_add(1, std::memory_order_relaxed);
				const int done = doneCount.fetch_add(1, std::memory_order_relaxed) + 1;
				if (progressCb) {
#pragma omp critical(covmc_param_progress)
					{
						if (!progressCb(done, N))
							cancelFlag.store(true);
					}
				}
				continue;
			}

			numPropagated.fetch_add(1, std::memory_order_relaxed);
			cloneConverged[static_cast<std::size_t>(k)] = 1;

#pragma omp critical(covmc_param_state)
			{
				cloneStates.push_back({R0.x, R0.y, R0.z, V0.x, V0.y, V0.z});
			}

			bool anyHit = false;
			for (std::size_t e = 0; e < Nev; ++e) {
				auto& sample = samples[e];
				const std::size_t slot = static_cast<std::size_t>(k) * Nev + e;
				missFlat[slot] = sample.missKm;
				if (sample.impact) {
					totalImpacts.fetch_add(1, std::memory_order_relaxed);
					anyHit = true;
					impactFlat[slot] = 1;
					sample.impact->cloneIndex = k;
					sample.impact->eventIndex = static_cast<int>(e);
					corridorByClone[static_cast<std::size_t>(k)].push_back(std::move(*sample.impact));
				}
			}

			if (anyHit)
				numAnyImpact.fetch_add(1, std::memory_order_relaxed);

			const int done = doneCount.fetch_add(1, std::memory_order_relaxed) + 1;
			if (progressCb) {
#pragma omp critical(covmc_param_progress)
				{
					if (!progressCb(done, N))
						cancelFlag.store(true);
				}
			}
		}

		for (int k = 0; k < N; ++k) {
			if (!cloneConverged[static_cast<std::size_t>(k)])
				continue;
			for (std::size_t e = 0; e < Nev; ++e) {
				const std::size_t slot = static_cast<std::size_t>(k) * Nev + e;
				const double miss = missFlat[slot];
				if (std::isnan(miss))
					continue;
				missAccum[e].push_back(miss);
				if (miss < out.minMissPerEventKm[e])
					out.minMissPerEventKm[e] = miss;
				if (impactFlat[slot])
					++out.hitsPerEvent[e];
			}
			for (ImpactCloneInfo& info : corridorByClone[static_cast<std::size_t>(k)]) {
				out.corridor.push_back(std::move(info));
				totalCorridorEntries.fetch_add(1, std::memory_order_relaxed);
			}
		}

		out.numConverged = numPropagated.load();
		out.numFailedFit = numFailed.load();
		out.numAnyImpactClones = numAnyImpact.load();
		out.cancelled = cancelFlag.load();

		if (out.numConverged > 0) {
			const double n = static_cast<double>(out.numConverged);
			const double pAny = static_cast<double>(out.numAnyImpactClones) / n;
			out.aggregateImpactProb = pAny;
			const double binom = std::max(pAny * (1.0 - pAny), 3.0 / n);
			out.aggregateImpactSigma = std::sqrt(binom / n);
		}

		for (std::size_t e = 0; e < Nev; ++e) {
			auto& v = missAccum[e];
			if (v.empty())
				continue;

			std::sort(v.begin(), v.end());
			out.medianMissPerEventKm[e] = v[v.size() / 2];
		}

		for (auto& m : out.minMissPerEventKm) {
			if (m == std::numeric_limits<double>::max())
				m = std::numeric_limits<double>::quiet_NaN();
		}

		const Vec3D centerR(parameterCenter(0), parameterCenter(1), parameterCenter(2));
		const Vec3D centerV(parameterCenter(3), parameterCenter(4), parameterCenter(5));
		const NonGravModel centerNg = nonGravFromPhysicalParameters(nominalNonGrav, parameterCenter);
		assignSurfacePoint(out, centerR, centerV, nominalOrbit.getEpoch(), &centerNg, propagator, events, eph,
		                   cfg.altitudeKm);

		const double convPct = (N > 0) ? 100.0 * static_cast<double>(out.numConverged) / static_cast<double>(N) : 0.0;

		LOG_INFO << "[ObsMC] N=" << N << "  converged=" << out.numConverged << " (" << std::fixed
		         << std::setprecision(0) << convPct << "%)"
		         << "  any-impact=" << out.numAnyImpactClones << "  P(any)=" << std::scientific << std::setprecision(3)
		         << out.aggregateImpactProb << " +/- " << out.aggregateImpactSigma;

		for (std::size_t e = 0; e < Nev; ++e) {
			LOG_INFO << "[ObsMC]   event[" << e << "] minMiss=" << std::fixed << std::setprecision(0)
			         << out.minMissPerEventKm[e] << " km"
			         << "  medianMiss=" << out.medianMissPerEventKm[e] << " km"
			         << "  hits=" << out.hitsPerEvent[e] << "/" << out.numConverged;
		}

		LOG_INFO << "[ObsMC] DIAG-XING summary: total impacts=" << totalImpacts.load()
		         << "  corridor entries added=" << totalCorridorEntries.load();

		logCorridorStatistics(out, events);

		if (cloneStates.size() >= 2) {
			std::array<double, 6> mean{};
			for (const auto& s : cloneStates) {
				for (int i = 0; i < 6; ++i)
					mean[i] += s[i];
			}
			for (int i = 0; i < 6; ++i)
				mean[i] /= static_cast<double>(cloneStates.size());

			std::array<double, 6> var{};
			for (const auto& s : cloneStates) {
				for (int i = 0; i < 6; ++i) {
					const double d = s[i] - mean[i];
					var[i] += d * d;
				}
			}

			const double denom = static_cast<double>(cloneStates.size() - 1);
			std::ostringstream sigStr;
			sigStr << std::scientific << std::setprecision(3);
			for (int i = 0; i < 6; ++i)
				sigStr << " " << std::sqrt(var[i] / denom);
			LOG_INFO << "[ObsMC]   empirical clone sigma (R AU, V AU/d):" << sigStr.str();
		}

		return out;
	}

	ObsMcResult runParameterCovarianceMonteCarlo(const Orbit& nominalOrbit,
	                                             const FitParameterCovariance& parameterCovariance,
	                                             const NonGravModel& nominalNonGrav, const NBodyPropagator& propagator,
	                                             std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg,
	                                             IEphemerisSource* eph, std::function<bool(int, int)> progressCb)
	{
		const VecX center = nominalPhysicalParameters(nominalOrbit, nominalNonGrav);
		return runParameterCovarianceMonteCarloAtCenter(nominalOrbit, center, parameterCovariance, nominalNonGrav,
		                                                propagator, events, cfg, eph, std::move(progressCb));
	}

	ObsMcResult runObservationMonteCarlo(const Orbit& nominalOrbit, std::span<const ObsNode> observations,
	                                     const FitConfig& fitCfg, const NBodyPropagator& propagator,
	                                     std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg,
	                                     IEphemerisSource* eph, std::function<bool(int, int)> progressCb)
	{
		if (!eph)
			eph = getDefaultEphemerisSource();
		validateAltitude(cfg.altitudeKm);
		ObsMcResult out;
		out.earthFrameModel = earthBodyFixedFrameModel();
		out.altitudeKm = cfg.altitudeKm;
		const int N = std::max(0, cfg.numClones);
		const std::size_t Nev = events.size();

		out.numAttempted = N;
		out.hitsPerEvent.assign(Nev, 0);
		out.minMissPerEventKm.assign(Nev, std::numeric_limits<double>::max());
		out.medianMissPerEventKm.assign(Nev, 0.0);

		if (N == 0 || observations.empty() || events.empty())
			return out;

		{
			const Vec3D R = nominalOrbit.getR();
			const Vec3D V = nominalOrbit.getV();
			LOG_INFO << "[ObsMC] start  rngSeed=" << cfg.rngSeed << "  N=" << N << "  sigmaMult=" << cfg.noiseSigmaMult
			         << "  Nobs=" << observations.size() << "  Nev=" << events.size() << "  epoch=" << std::fixed
			         << std::setprecision(9) << nominalOrbit.getEpoch() << "  R=(" << R.x << "," << R.y << "," << R.z
			         << ")"
			         << "  V=(" << V.x << "," << V.y << "," << V.z << ")";
		}

		const std::size_t Nobs = observations.size();
		std::vector<double> sigmaRaRad(Nobs), sigmaDecRad(Nobs);
		constexpr double kDefaultSigmaRad = Constants::ARCSEC_RAD;
		for (std::size_t i = 0; i < Nobs; ++i) {
			double sRa = kDefaultSigmaRad;
			double sDec = kDefaultSigmaRad;
			const auto uncertainty = obs::lookupUncertainty(observations[i], observations[i].JD, fitCfg.weightProfile);
			sRa = uncertainty.sigmaRaRad;
			sDec = uncertainty.sigmaDecRad;
			if (obs::usesObserverSigma(observations[i], fitCfg.obsWeighting.useObserverSigmas)) {
				sRa = observations[i].rmsRaArcsec * Constants::ARCSEC_RAD;
				sDec = observations[i].rmsDecArcsec * Constants::ARCSEC_RAD;
			}
			sigmaRaRad[i] = sRa * cfg.noiseSigmaMult;
			sigmaDecRad[i] = sDec * cfg.noiseSigmaMult;
		}

		std::vector<double> eventGamma(Nev);
		for (std::size_t e = 0; e < Nev; ++e) {
			eventGamma[e] = gravitationalFocusing(events[e].vInfKmS, events[e].planet);
			const double rPl = planetRadiusKm(events[e].planet);
			LOG_INFO << "[ObsMC] event[" << e << "] planet=" << static_cast<int>(events[e].planet)
			         << "  v_inf=" << std::fixed << std::setprecision(2) << events[e].vInfKmS << " km/s"
			         << "  γ=" << std::setprecision(4) << eventGamma[e] << "  R=" << std::setprecision(0) << rPl
			         << " km"
			         << "  surface-radius=" << rPl << " km"
			         << "  (γR=" << eventGamma[e] * rPl << " km, info only)";
		}

		std::vector<std::vector<double>> missAccum(Nev);
		std::vector<std::vector<double>> qAccum(Nev);

		std::vector<std::array<double, 6>> cloneStates;
		std::vector<int> cloneIterations;
		std::vector<double> cloneCosts;

		std::atomic<int> numConverged{0};
		std::atomic<int> numFailedFit{0};
		std::atomic<int> numAnyImpact{0};
		std::atomic<int> doneCount{0};
		std::atomic<bool> cancelFlag{false};

		const Orbit nominalForClone(nominalOrbit);

		const FitRunOptions cloneRunOptions{
		    .selection = FitSelectionMode::Fixed,
		    .uncertainty = FitUncertaintyMode::Skip,
		};

#pragma omp parallel for schedule(dynamic, 4)
		for (int k = 0; k < N; ++k) {
			if (cancelFlag.load(std::memory_order_relaxed))
				continue;

			const std::uint64_t seed = cfg.rngSeed ^ (0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(k + 1));
			std::mt19937_64 rng(seed);
			std::normal_distribution<double> N01(0.0, 1.0);

			std::vector<ObsNode> noisy(observations.begin(), observations.end());
			for (std::size_t i = 0; i < Nobs; ++i) {
				if (noisy[i].radar) {
					auto& radar = *noisy[i].radar;
					if (radar.delaySeconds)
						radar.delaySeconds->value += N01(rng) * radar.delaySeconds->sigma * cfg.noiseSigmaMult;
					if (radar.dopplerHz)
						radar.dopplerHz->value += N01(rng) * radar.dopplerHz->sigma * cfg.noiseSigmaMult;
					continue;
				}
				const double cosD = std::max(std::abs(std::cos(noisy[i].DEC)), 1.0e-9);
				const double dra = N01(rng) * sigmaRaRad[i] / cosD;
				const double ddec = N01(rng) * sigmaDecRad[i];
				noisy[i].RA += dra;
				noisy[i].DEC += ddec;
				recomputeDirectionVector(noisy[i]);
			}

			ScopedLogSilence cloneSilence(LogLevel::Warning);

			auto cloneOrbit = std::make_shared<Orbit>(nominalForClone);
			FitSession fs(fitCfg);
			FitResult fr;
			try {
				fr = fs.optimize({.initialOrbit = cloneOrbit, .observations = noisy, .options = cloneRunOptions});
			}

			catch (...) {
				numFailedFit.fetch_add(1, std::memory_order_relaxed);
				continue;
			}

			if (!fr.converged || !fr.orbit) {
				numFailedFit.fetch_add(1, std::memory_order_relaxed);
				continue;
			}

			const Vec3D R0 = fr.orbit->getR();
			const Vec3D V0 = fr.orbit->getV() * kGaussK;
			const double epoch = fr.orbit->getEpoch();

			std::vector<EncounterSample> samples;
			try {
				samples = sampleEncounters({R0, V0, epoch, fr.finalNonGrav}, events, propagator, eph, cfg.altitudeKm);
			} catch (...) {
				numFailedFit.fetch_add(1, std::memory_order_relaxed);
				continue;
			}
			numConverged.fetch_add(1, std::memory_order_relaxed);

#pragma omp critical(obsmc_state)
			{
				cloneStates.push_back({R0.x, R0.y, R0.z, V0.x, V0.y, V0.z});
				cloneIterations.push_back(fr.iterations);
				cloneCosts.push_back(fr.cost);
			}

			bool anyHit = false;
			for (std::size_t e = 0; e < Nev; ++e) {
				auto& sample = samples[e];
#pragma omp critical(obsmc_accum)
				{
					if (std::isfinite(sample.missKm)) {
						missAccum[e].push_back(sample.missKm);
						qAccum[e].push_back(sample.conicPericenterKm);
						out.minMissPerEventKm[e] = std::min(out.minMissPerEventKm[e], sample.missKm);
					}
					if (sample.impact) {
						++out.hitsPerEvent[e];
						anyHit = true;
						sample.impact->cloneIndex = k;
						sample.impact->eventIndex = static_cast<int>(e);
						out.corridor.push_back(std::move(*sample.impact));
					}
				}
			}
			if (anyHit)
				numAnyImpact.fetch_add(1, std::memory_order_relaxed);

			const int done = doneCount.fetch_add(1, std::memory_order_relaxed) + 1;
			if (progressCb) {
#pragma omp critical(obsmc_progress)
				{
					if (!progressCb(done, N))
						cancelFlag.store(true);
				}
			}
		}

		out.numConverged = numConverged.load();
		out.numFailedFit = numFailedFit.load();
		out.numAnyImpactClones = numAnyImpact.load();
		out.cancelled = cancelFlag.load();

		if (out.numConverged > 0) {
			const double n = static_cast<double>(out.numConverged);
			const double pAny = static_cast<double>(out.numAnyImpactClones) / n;
			out.aggregateImpactProb = pAny;
			const double q = std::max(pAny * (1.0 - pAny), 3.0 / n);
			out.aggregateImpactSigma = std::sqrt(q / n);
		}

		for (std::size_t e = 0; e < Nev; ++e) {
			auto& v = missAccum[e];
			if (v.empty())
				continue;

			std::sort(v.begin(), v.end());
			out.medianMissPerEventKm[e] = v[v.size() / 2];
		}

		for (auto& m : out.minMissPerEventKm)
			if (m == std::numeric_limits<double>::max())
				m = std::numeric_limits<double>::quiet_NaN();

		const double convPct = (N > 0) ? 100.0 * static_cast<double>(out.numConverged) / static_cast<double>(N) : 0.0;

		LOG_INFO << "[ObsMC] N=" << N << "  converged=" << out.numConverged << " (" << std::fixed
		         << std::setprecision(0) << convPct << "%)"
		         << "  any-impact=" << out.numAnyImpactClones << "  P(any)=" << std::scientific << std::setprecision(3)
		         << out.aggregateImpactProb << " ± " << out.aggregateImpactSigma;

		for (std::size_t e = 0; e < Nev; ++e) {
			LOG_INFO << "[ObsMC]   event[" << e << "] minMiss=" << std::fixed << std::setprecision(0)
			         << out.minMissPerEventKm[e] << " km"
			         << "  medianMiss=" << out.medianMissPerEventKm[e] << " km"
			         << "  hits=" << out.hitsPerEvent[e] << "/" << out.numConverged;
		}

		logCorridorStatistics(out, events);

		for (std::size_t e = 0; e < Nev; ++e) {
			auto& v = qAccum[e];
			if (v.empty())
				continue;
			const double R = planetRadiusKm(events[e].planet);
			const double gR = eventGamma[e] * R;

			int below_R = 0;
			int in_grazing = 0;
			int in_2R = 0;
			int in_5R = 0;
			int in_10R = 0;
			int in_100R = 0;
			int beyond_100R = 0;

			for (double q : v) {
				if (q < R)
					++below_R;

				else if (q < gR)
					++in_grazing;

				else if (q < 2 * R)
					++in_2R;

				else if (q < 5 * R)
					++in_5R;

				else if (q < 10 * R)
					++in_10R;

				else if (q < 100 * R)
					++in_100R;

				else
					++beyond_100R;
			}

			std::sort(v.begin(), v.end());
			const std::size_t n = v.size();
			auto pct = [&](double p) -> double {
				std::size_t k = static_cast<std::size_t>(p * (n - 1));
				return v[k];
			};

			LOG_INFO << "[ObsMC]   event[" << e << "] q-distribution"
			         << " (R=" << std::fixed << std::setprecision(0) << R << " km, γR=" << gR << " km)";

			LOG_INFO << "[ObsMC]     q < R         : " << below_R << "    R   ≤ q < γR    : " << in_grazing
			         << "    γR  ≤ q < 2R    : " << in_2R;

			LOG_INFO << "[ObsMC]     2R  ≤ q < 5R  : " << in_5R << "    5R  ≤ q < 10R   : " << in_10R
			         << "    10R ≤ q < 100R  : " << in_100R << "    q ≥ 100R        : " << beyond_100R;

			LOG_INFO << "[ObsMC]    q percentiles:"
			         << " p1=" << std::fixed << std::setprecision(0) << pct(0.01) << " km"
			         << " p10=" << pct(0.10) << " km"
			         << " p50=" << pct(0.50) << " km"
			         << " p90=" << pct(0.90) << " km"
			         << " p99=" << pct(0.99) << " km";
		}

		for (std::size_t e = 0; e < Nev; ++e) {
			int nValid = 0;
			int nInvalid = 0;
			double minJd = std::numeric_limits<double>::max();
			double maxJd = std::numeric_limits<double>::lowest();

			int nAll = 0;

			double minLat = std::numeric_limits<double>::max();
			double maxLat = std::numeric_limits<double>::lowest();
			double minLon = std::numeric_limits<double>::max();
			double maxLon = std::numeric_limits<double>::lowest();

			constexpr int kLonBins = 24;
			constexpr int kLatBins = 12;
			constexpr double kLonBinDeg = 360.0 / kLonBins;
			constexpr double kLatBinDeg = 180.0 / kLatBins;
			std::array<int, kLonBins> lonHist{};
			std::array<int, kLatBins> latHist{};

			for (const auto& info : out.corridor) {
				if (info.eventIndex != static_cast<int>(e))
					continue;

				++nAll;
				if (!info.bodyFixedValid) {
					++nInvalid;
					continue;
				}

				++nValid;
				minJd = std::min(minJd, info.tcaJd);
				maxJd = std::max(maxJd, info.tcaJd);
				minLat = std::min(minLat, info.latitudeDeg);
				maxLat = std::max(maxLat, info.latitudeDeg);
				minLon = std::min(minLon, info.longitudeDeg);
				maxLon = std::max(maxLon, info.longitudeDeg);

				int lonIdx = static_cast<int>(std::floor((info.longitudeDeg + 180.0) / kLonBinDeg));
				lonIdx = std::clamp(lonIdx, 0, kLonBins - 1);
				++lonHist[lonIdx];

				int latIdx = static_cast<int>(std::floor((info.latitudeDeg + 90.0) / kLatBinDeg));
				latIdx = std::clamp(latIdx, 0, kLatBins - 1);
				++latHist[latIdx];
			}

			if (nValid == 0) {
				LOG_INFO << "[ObsMC]   event[" << e << "] corridor: no valid lat/lon entries"
				         << "  (total=" << nAll << "  invalid=" << nInvalid << ")";
				continue;
			}

			const double tcaSpreadHours = (maxJd - minJd) * 24.0;

			LOG_INFO << "[ObsMC]   event[" << e << "] corridor counts: "
			         << "total=" << nAll << "  bodyFixedValid=" << nValid << "  invalid=" << nInvalid << "  TCA spread "
			         << std::fixed << std::setprecision(2) << tcaSpreadHours << " h";

			LOG_INFO << "[ObsMC]   event[" << e << "] lat range " << std::fixed << std::setprecision(3) << "[" << minLat
			         << ", " << maxLat << "]°  "
			         << "lon range [" << minLon << ", " << maxLon << "]°";

			{
				std::ostringstream h;
				h << "[ObsMC]   event[" << e << "] lon hist (15°): ";
				for (int b = 0; b < kLonBins; ++b) {
					if (lonHist[b] == 0)
						continue;
					const double lo = -180.0 + b * kLonBinDeg;
					h << "[" << std::fixed << std::setprecision(0) << lo << ".." << (lo + kLonBinDeg)
					  << "]=" << lonHist[b] << " ";
				}
				LOG_INFO << h.str();
			}

			{
				std::ostringstream h;
				h << "[ObsMC]   event[" << e << "] lat hist (15°): ";
				for (int b = 0; b < kLatBins; ++b) {
					if (latHist[b] == 0)
						continue;
					const double lo = -90.0 + b * kLatBinDeg;
					h << "[" << std::fixed << std::setprecision(0) << lo << ".." << (lo + kLatBinDeg)
					  << "]=" << latHist[b] << " ";
				}
				LOG_INFO << h.str();
			}
		}

		assignNominalSurfacePoint(out, nominalOrbit, propagator, events, eph, cfg.altitudeKm);

		if (cloneStates.size() >= 2) {
			std::array<double, 6> mean{};
			for (const auto& s : cloneStates)
				for (int i = 0; i < 6; ++i)
					mean[i] += s[i];

			for (int i = 0; i < 6; ++i)
				mean[i] /= static_cast<double>(cloneStates.size());

			std::array<double, 6> var{};
			for (const auto& s : cloneStates)
				for (int i = 0; i < 6; ++i) {
					const double d = s[i] - mean[i];
					var[i] += d * d;
				}

			const double denom = static_cast<double>(cloneStates.size() - 1);
			std::array<double, 6> sig{};
			for (int i = 0; i < 6; ++i)
				sig[i] = std::sqrt(var[i] / denom);

			std::ostringstream sigStr;
			sigStr << std::scientific << std::setprecision(3);
			for (int i = 0; i < 6; ++i)
				sigStr << " " << sig[i];

			LOG_INFO << "[ObsMC]   empirical clone σ (R AU, V AU/d):" << sigStr.str();

			std::ostringstream meanStr;
			meanStr << std::scientific << std::setprecision(3);

			for (int i = 0; i < 6; ++i)
				meanStr << " " << mean[i];

			LOG_INFO << "[ObsMC]   empirical clone mean state:" << meanStr.str();
		}

		if (!cloneIterations.empty()) {
			auto iterMin = *std::min_element(cloneIterations.begin(), cloneIterations.end());
			auto iterMax = *std::max_element(cloneIterations.begin(), cloneIterations.end());
			double iterMean = 0.0;

			for (int it : cloneIterations)
				iterMean += static_cast<double>(it);

			iterMean /= static_cast<double>(cloneIterations.size());
			double costMin = *std::min_element(cloneCosts.begin(), cloneCosts.end());
			double costMax = *std::max_element(cloneCosts.begin(), cloneCosts.end());
			double costMean = 0.0;
			for (double c : cloneCosts)
				costMean += c;

			costMean /= static_cast<double>(cloneCosts.size());
			LOG_INFO << "[ObsMC]   per-clone iters: min=" << iterMin << "  mean=" << std::fixed << std::setprecision(1)
			         << iterMean << "  max=" << iterMax << "  cost: min=" << std::scientific << std::setprecision(3)
			         << costMin << "  mean=" << costMean << "  max=" << costMax;
		}

		if (!out.isReliable()) {
			LOG_INFO << "[ObsMC]   WARNING: convergence < 50% — arc too short, "
			            "results unreliable.";
		}

		return out;
	}
}  // namespace orbital_core::hazard
