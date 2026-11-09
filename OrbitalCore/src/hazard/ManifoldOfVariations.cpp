#include "orbital_core/hazard/ManifoldOfVariations.h"
#include "mov/MovCoordinates.h"
#include "fit/FitObservations.h"
#include "fit/OrbitFitModel.h"
#include "fit/OrbitFitProblem.h"
#include "obs/ObserverGeometry.h"
#include "hazard/McBoundaryScan.h"
#include "hazard/BoundaryScan.h"
#include "orbital_core/optimizers/EigenLevenbergMarquardt.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "logger/Logger.h"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <mutex>
#include <numeric>

namespace orbital_core::hazard {
	namespace {
		using mov::Matrix6;
		using mov::Vector6;

		Vec3D vec(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

		std::array<double, 6> array(const Vector6& v)
		{
			std::array<double, 6> result;
			std::copy(v.data(), v.data() + 6, result.begin());
			return result;
		}

		class ConstrainedFit final : public optimizers::ILSQProblem {
		   public:
			ConstrainedFit(const mov::Chart& chart, const Vector6& initial,
			               const std::vector<obs::PreparedObservation>& observations,
			               const OrbitOptimizerConfig& config, double slopeG)
			    : chart_(chart),
			      fixed_(initial),
			      model_(orbit(chart, initial), observations, config),
			      problem_(model_),
			      slopeG_(slopeG)
			{
				problem_.beginPass(6);
			}

			static Orbit_p orbit(const mov::Chart& chart, const Vector6& att)
			{
				const auto state = chart.state(att);
				auto result = std::make_shared<Orbit>();
				result->setEpochMjd(chart.epoch);
				result->setR(vec(state.head<3>()));
				result->setV(vec(state.tail<3>() / Constants::GAUSS_K));
				return result;
			}

			int numParams() const noexcept override { return 4; }

			int numResiduals() const noexcept override { return model_.numResiduals(); }

			Vector6 attributable(const std::vector<double>& params) const
			{
				auto result = fixed_;
				for (int i = 0; i < 4; ++i)
					result(i) = params[i];
				return result;
			}

			std::vector<double> stateParameters(const std::vector<double>& params) const
			{
				Vector6 state = chart_.state(attributable(params));
				state.tail<3>() /= Constants::GAUSS_K;
				return {state.data(), state.data() + 6};
			}

			void residuals(const std::vector<double>& params, std::vector<double>& out) const override
			{
				problem_.residuals(stateParameters(params), out);
			}

			Eigen::MatrixXd design(const std::vector<double>& params) const
			{
				Eigen::MatrixXd j;
				problem_.jacobian(stateParameters(params), j);
				Matrix6 d = chart_.derivative(attributable(params));
				d.bottomRows<3>() /= Constants::GAUSS_K;
				return j * d;
			}

			void jacobian(const std::vector<double>& params, Eigen::MatrixXd& out) const override
			{
				out = design(params).leftCols(4);
			}

			bool acceptParams(const std::vector<double>& params) const override
			{
				return params.size() == 4 &&
				       std::all_of(params.begin(), params.end(), [](double x) { return std::isfinite(x); }) &&
				       std::abs(params[1]) < Constants::PI / 2 - 1e-8 && std::abs(params[2]) < 10 &&
				       std::abs(params[3]) < 10;
			}

			double displayCost(const std::vector<double>& params) const override
			{
				std::vector<double> r;
				residuals(params, r);
				return std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
			}

			std::optional<double> magnitude(const std::vector<double>& params) const
			{
				const auto states = model_.computeLightTimePositions(stateParameters(params), NonGravModel{});
				const auto& obs = model_.observations();
				if (states.size() != obs.size())
					throw std::runtime_error("MOV photometric geometry unavailable");
				std::vector<Residual> photometryRows;
				for (std::size_t i = 0; i < obs.size(); ++i) {
					const auto relative = states[i] - obs[i].earthR - obs[i].observerOffsetEclAu;
					Residual row;
					row.mag = obs[i].mag;
					row.magBand = obs[i].magBand;
					row.magSigma = obs[i].magSigma;
					row.helioDistAu = states[i].magnitude();
					row.geoDistAu = relative.magnitude();
					row.phaseAngleRad =
					    std::acos(std::clamp(states[i].dot(relative) / (row.helioDistAu * row.geoDistAu), -1.0, 1.0));
					photometryRows.push_back(row);
				}
				return photometry::fit(photometryRows, slopeG_).absoluteMagnitude;
			}

		   private:
			const mov::Chart& chart_;
			Vector6 fixed_;
			fit::OrbitFitModel model_;
			mutable fit::OrbitFitProblem problem_;
			double slopeG_;
		};

		struct GridBounds {
			std::vector<std::pair<double, double>> ranges;
			double rateMin{}, rateMax{};
			bool logarithmic{true};
		};

		std::vector<MovNode> grid(const GridBounds& bounds, int ranges, int rates)
		{
			std::vector<MovNode> nodes;
			for (std::size_t component = 0; component < bounds.ranges.size(); ++component) {
				const auto [low, high] = bounds.ranges[component];
				const double first = bounds.logarithmic ? std::log10(low) : low;
				const double last = bounds.logarithmic ? std::log10(high) : high;
				const double area = (last - first) / (ranges - 1) * (bounds.rateMax - bounds.rateMin) / (rates - 1);
				for (int i = 0; i < ranges; ++i)
					for (int j = 0; j < rates; ++j) {
						MovNode n;
						n.row = static_cast<int>(component) * ranges + i;
						n.column = j;
						const double coordinate = first + (last - first) * i / (ranges - 1);
						n.attributable[4] = bounds.logarithmic ? std::pow(10.0, coordinate) : coordinate;
						n.attributable[5] = bounds.rateMin + (bounds.rateMax - bounds.rateMin) * j / (rates - 1);
						n.jacobianSampling = area * (bounds.logarithmic ? std::log(10.0) * n.attributable[4] : 1);
						n.outerBoundary = i == 0 || i + 1 == ranges || j == 0 || j + 1 == rates;
						nodes.push_back(n);
					}
			}
			return nodes;
		}
	}  // namespace

	MovResult computeManifoldOfVariations(std::span<const ObsNode> observations, const FitConfig& settings,
	                                      const MovConfig& config, const FitResult* nominal, std::stop_token cancel,
	                                      std::function<void(int, int)> progress)
	{
		if (!settings.nbody.enabled || !settings.useWeights || settings.nonGrav.enabled)
			throw std::invalid_argument("MOV requires weighted optical observations and gravity-only N-body dynamics");
		if ((config.sampling != MovSampling::Automatic && config.sampling != MovSampling::Grid &&
		     config.sampling != MovSampling::Cobweb) ||
		    config.preliminaryNodes < 8 || config.preliminaryNodes > 250 || config.rangeNodes < 8 ||
		    config.rangeNodes > 250 || config.rateNodes < 8 || config.rateNodes > 250 || config.workers < 1 ||
		    config.workers > 1024 || !std::isfinite(config.startMjdTdb) || !std::isfinite(config.endMjdTdb) ||
		    config.endMjdTdb <= config.startMjdTdb || !std::isfinite(config.sigmaMax) || config.sigmaMax <= 0 ||
		    config.sigmaMax > 10 || !std::isfinite(config.absoluteMagnitudeMax) || !std::isfinite(config.altitudeKm) ||
		    config.altitudeKm < 0 ||
		    (config.target != Planet::Earth && config.target != Planet::Moon && config.target != Planet::Mars))
			throw std::invalid_argument("Invalid MOV configuration");
		auto* eph = getDefaultEphemerisSource();
		if (!eph || !eph->isAvailable())
			throw std::runtime_error("Planetary ephemerides unavailable");
		if (cancel.stop_requested()) {
			MovResult r;
			r.config = config;
			r.cancelled = true;
			return r;
		}
		OrbitOptimizerConfig fitConfig;
		fitConfig.nbody = settings.nbody;
		fitConfig.useWeights = true;
		fitConfig.weightProfile = settings.weightProfile;
		fitConfig.obsWeighting = settings.obsWeighting;
		fitConfig.logDiagnostics = false;
		fitConfig.logCacheStats = false;
		auto prepared = fit::prepareFitObservations(observations, fitConfig);
		const auto attributable = mov::fitAttributable(prepared);
		MovResult result;
		result.config = config;
		result.epochMjdTdb = attributable.epoch;
		result.curvatureSnr = attributable.curvatureSnr;
		std::copy(attributable.angles.data(), attributable.angles.data() + 4, result.attributable.begin());
		result.observerCode = attributable.observer.obsCode;
		if (config.startMjdTdb < result.epochMjdTdb || config.endMjdTdb - result.epochMjdTdb > 30)
			throw std::invalid_argument(
			    "MOV interval must start at/after its attributable epoch and end within 30 days");
		mov::Chart chart{attributable.epoch, [eph, observer = attributable.observer](double t) {
			                 const auto earth = eph->queryMjd(EphemerisPlanet::Earth, t);
			                 const auto offset =
			                     obs::groundObserverStateEclAuMjd(t, time::fromMjdTdb(t, time::Scale::UTC), observer);
			                 return NBodyPropagator::StateVel{earth.R + offset.offsetEclAu,
			                                                  earth.V + offset.velocityEclAuPerDay};
		                 }};
		const auto region = mov::admissibleRegion(chart, attributable.angles);
		NBodyPropagator propagator(IntegratorKind::RA15);
		propagator.setStepSize(settings.nbody.stepSize);
		propagator.setPlanetMask(settings.nbody.planetMask);
		propagator.setUseSb441(settings.nbody.useSb441);
		propagator.setTargetNaifId(settings.nbody.targetNaifId);
		propagator.setInactiveInteriorMassPolicy(settings.nbody.inactiveInteriorMassPolicy);
		propagator = propagator.forEncounterMonitoring(*eph);
		const auto frame = bodyFixedTransformMjd(config.target, config.startMjdTdb, eph);
		if (!frame.valid)
			throw std::runtime_error("MOV target orientation unavailable");
		result.bodyFixedFrame = frame.frame;
		result.bodyFixedModel = config.target == Planet::Earth ? earthBodyFixedFrameModel() : std::string(frame.model);
		Vector6 initial;
		initial << attributable.angles, 1, 0;
		Eigen::Matrix2d rangeCovariance;
		std::string nominalUnavailable;
		if (!nominal || !nominal->fitEpochOrbit)
			nominalUnavailable = "no nominal orbit supplied for the selected observations";
		else if (!nominal->converged || nominal->cancelled)
			nominalUnavailable = "the nominal fit for the selected observations has not converged";
		else if (nominal->uncertainty.fitParameterCovariance.rows() != 6 ||
		         nominal->uncertainty.fitParameterCovariance.cols() != 6)
			nominalUnavailable = "the nominal fit has no full 6 x 6 covariance";
		else if (!(attributable.curvatureSnr > 3))
			nominalUnavailable = "curvature SNR = " + std::to_string(attributable.curvatureSnr) + " (required > 3)";
		const bool useNominal = config.sampling != MovSampling::Grid && nominalUnavailable.empty();
		if (config.sampling == MovSampling::Cobweb && !useNominal)
			throw std::invalid_argument("Cobweb unavailable: " + nominalUnavailable +
			                            ". Use Automatic or admissible-region grid.");
		if (useNominal) {
			const auto& orbit = *nominal->fitEpochOrbit;
			auto p = propagator;
			const auto states =
			    p.propagateAllWithStmAndSensMjd(orbit.getR(), orbit.getV() * Constants::GAUSS_K, orbit.epochMjd(),
			                                    {result.epochMjdTdb}, NonGravModel{}, nullptr, cancel);
			if (states.size() != 1) {
				if (cancel.stop_requested()) {
					result.cancelled = true;
					return result;
				}
				throw std::runtime_error("MOV nominal covariance transport did not complete");
			}
			Vector6 x;
			x << states[0].R.x, states[0].R.y, states[0].R.z, states[0].V.x, states[0].V.y, states[0].V.z;
			initial = chart.inverse(x);
			if (!nominal->uncertainty.fitParameterCovariance.allFinite() ||
			    !nominal->uncertainty.fitParameterCovariance.isApprox(
			        nominal->uncertainty.fitParameterCovariance.transpose(), 1e-10))
				throw std::invalid_argument("Invalid MOV nominal covariance");
			const Matrix6 inverse = chart.derivative(initial).inverse();
			const Matrix6 transport = inverse * states[0].Phi;
			const Matrix6 covariance = transport * nominal->uncertainty.fitParameterCovariance * transport.transpose();
			rangeCovariance = covariance.bottomRightCorner<2, 2>();
			const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(rangeCovariance);
			if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() <= 0)
				throw std::runtime_error("MOV range covariance is not positive definite");
		}
		result.sampling = useNominal ? MovSampling::Cobweb : MovSampling::Grid;
		std::mutex mutex;
		int done = 0, total = 0;
		std::exception_ptr callbackError;
		const auto fitNodes = [&](std::vector<MovNode>& nodes) {
			total += static_cast<int>(nodes.size());
			std::vector<std::pair<int, int>> rows;
			for (int first = 0; first < static_cast<int>(nodes.size());) {
				int last = first + 1;
				while (last < static_cast<int>(nodes.size()) && nodes[last].row == nodes[first].row)
					++last;
				rows.emplace_back(first, last);
				first = last;
			}
#pragma omp parallel for schedule(dynamic, 1) num_threads(config.workers)
			for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
				Eigen::Vector4d guess = initial.head<4>();
				for (int i = rows[row].first; i < rows[row].second; ++i) {
					auto& node = nodes[i];
					if (cancel.stop_requested()) {
						node.status = MovNodeStatus::Cancelled;
						continue;
					}
					try {
						if (!region.contains(node.attributable[4], node.attributable[5]))
							node.status = MovNodeStatus::OutsideRegion;
						else {
							Vector6 att = initial;
							att.head<4>() = guess;
							att(4) = node.attributable[4];
							att(5) = node.attributable[5];
							ScopedLogSilence silence(LogLevel::Warning);
							ConstrainedFit problem(chart, att, prepared, fitConfig, settings.photometry.asteroidG);
							const auto fit = optimizers::EigenLevenbergMarquardt{}.optimize(
							    {att.data(), att.data() + 4}, problem, [&] { return cancel.stop_requested(); });
							if (fit.cancelled)
								node.status = MovNodeStatus::Cancelled;
							else if (!fit.converged) {
								node.status = MovNodeStatus::FitFailed;
								node.error = "Constrained fit did not converge";
							} else {
								const Vector6 fitted = problem.attributable(fit.params);
								guess = fitted.head<4>();
								node.attributable = array(fitted);
								node.state = array(chart.state(fitted));
								node.chiSquare = problem.displayCost(fit.params);
								node.jacobianMov = mov::manifoldJacobian(problem.design(fit.params));
								node.absoluteMagnitude = problem.magnitude(fit.params);
								node.status = MovNodeStatus::Rejected;
								if (!node.absoluteMagnitude)
									node.error = "Absolute magnitude unavailable";
							}
						}
					} catch (const std::exception& e) {
						node.status = MovNodeStatus::FitFailed;
						node.error = e.what();
					}
					std::lock_guard lock(mutex);
					++done;
					if (progress && !callbackError)
						try {
							progress(done, total);
						} catch (...) {
							callbackError = std::current_exception();
						}
				}
			}
			if (callbackError)
				std::rethrow_exception(callbackError);
			double minimum = std::numeric_limits<double>::infinity();
			for (const auto& n : nodes)
				if (n.status == MovNodeStatus::Rejected && n.absoluteMagnitude &&
				    *n.absoluteMagnitude <= config.absoluteMagnitudeMax)
					minimum = std::min(minimum, n.chiSquare);
			for (auto& n : nodes)
				if (std::isfinite(n.chiSquare) && std::isfinite(minimum))
					n.chi = std::sqrt(std::max(0.0, n.chiSquare - minimum));
			return minimum;
		};
		const auto eligible = [&](const MovNode& n) {
			return n.status == MovNodeStatus::Rejected && std::isfinite(n.chi) && n.chi <= config.sigmaMax &&
			       n.absoluteMagnitude && *n.absoluteMagnitude <= config.absoluteMagnitudeMax;
		};
		if (useNominal) {
			MovNode centre;
			centre.row = -1;
			centre.attributable = array(initial);
			result.nodes.push_back(centre);
			const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(rangeCovariance);
			const Eigen::Matrix2d axes = eig.eigenvectors() * eig.eigenvalues().cwiseSqrt().asDiagonal();
			for (int i = 0; i < config.rangeNodes; ++i)
				for (int j = 1; j <= config.rateNodes; ++j) {
					const double radius = config.sigmaMax * j / config.rateNodes,
					             angle = 2 * Constants::PI * (i + 1) / config.rangeNodes;
					const Eigen::Vector2d polar =
					    initial.tail<2>() + radius * axes * Eigen::Vector2d(std::cos(angle), std::sin(angle));
					MovNode node;
					node.row = i;
					node.column = j;
					node.attributable[4] = polar.x();
					node.attributable[5] = polar.y();
					node.jacobianSampling = radius * std::abs(axes.determinant()) * config.sigmaMax / config.rateNodes *
					                        2 * Constants::PI / config.rangeNodes;
					node.outerBoundary = j == config.rateNodes;
					result.nodes.push_back(node);
				}
		} else {
			GridBounds bounds{region.ranges, region.rateMin, region.rateMax,
			                  region.ranges.front().second < std::sqrt(10.0)};
			result.preliminaryLogRange = bounds.logarithmic;
			result.preliminary = grid(bounds, config.preliminaryNodes, config.preliminaryNodes);
			fitNodes(result.preliminary);
			if (cancel.stop_requested()) {
				result.cancelled = true;
				return result;
			}
			std::vector<MovNode> accepted;
			for (const auto& n : result.preliminary)
				if (eligible(n))
					accepted.push_back(n);
			if (accepted.empty())
				throw std::runtime_error("No admissible converged preliminary MOV nodes");
			double nearEarthWeight = 0, totalWeight = 0;
			for (const auto& n : accepted) {
				const double weight = std::exp(-n.chi * n.chi / 2) * n.jacobianMov * n.jacobianSampling;
				const auto& x = n.state;
				const Vec3D r{x[0], x[1], x[2]}, v{x[3], x[4], x[5]};
				const Vec3D h = r.cross(v), eccentricity = v.cross(h) / Constants::GM_SUN - r / r.magnitude();
				const double perihelion = h.dot(h) / (Constants::GM_SUN * (1 + eccentricity.magnitude()));
				totalWeight += weight;
				if (perihelion < 1.3)
					nearEarthWeight += weight;
			}
			bounds.ranges.clear();
			bounds.rateMin = region.rateMax;
			bounds.rateMax = region.rateMin;
			const double rateStep = (region.rateMax - region.rateMin) / (config.preliminaryNodes - 1);
			for (auto interval : region.ranges) {
				double low = interval.second, high = interval.first;
				bool found = false;
				for (const auto& n : accepted)
					if (n.attributable[4] >= interval.first && n.attributable[4] <= interval.second) {
						low = std::min(low, n.attributable[4]);
						high = std::max(high, n.attributable[4]);
						found = true;
						bounds.rateMin = std::min(bounds.rateMin, n.attributable[5] - rateStep);
						bounds.rateMax = std::max(bounds.rateMax, n.attributable[5] + rateStep);
					}
				if (found) {
					if (bounds.logarithmic) {
						const double step =
						    std::pow(interval.second / interval.first, 1.0 / (config.preliminaryNodes - 1));
						bounds.ranges.emplace_back(std::max(interval.first, low / step),
						                           std::min(interval.second, high * step));
					} else {
						const double step = (interval.second - interval.first) / (config.preliminaryNodes - 1);
						bounds.ranges.emplace_back(std::max(interval.first, low - step),
						                           std::min(interval.second, high + step));
					}
				}
			}
			bounds.rateMin = std::max(bounds.rateMin, region.rateMin);
			bounds.rateMax = std::min(bounds.rateMax, region.rateMax);
			if (region.ranges.size() == 1)
				bounds.logarithmic = nearEarthWeight > totalWeight / 2;
			result.logRange = bounds.logarithmic;
			result.nodes = grid(bounds, config.rangeNodes, config.rateNodes);
		}
		fitNodes(result.nodes);
		if (cancel.stop_requested()) {
			result.cancelled = true;
			return result;
		}
		IntervalMcConfig scan;
		scan.target = config.target;
		scan.startMjdTdb = config.startMjdTdb;
		scan.endMjdTdb = config.endMjdTdb;
		scan.altitudeKm = config.altitudeKm;
		if (nominal && nominal->converged && !nominal->cancelled && nominal->fitEpochOrbit) {
			MovNominal trajectory;
			const auto& orbit = *nominal->fitEpochOrbit;
			const auto r = orbit.getR(), v = orbit.getV() * Constants::GAUSS_K;
			trajectory.epochMjdTdb = orbit.epochMjd();
			trajectory.state = {r.x, r.y, r.z, v.x, v.y, v.z};
			try {
				if (trajectory.epochMjdTdb > scan.startMjdTdb)
					throw std::runtime_error("Nominal epoch is later than the MOV scan start");
				auto p = propagator;
				const auto crossing = boundary::search(trajectory.state, trajectory.epochMjdTdb, nominal->finalNonGrav,
				                                       scan.target, scan.startMjdTdb, scan.endMjdTdb, scan.altitudeKm,
				                                       p, *eph, cancel);
				trajectory.entryBeforeInterval = crossing.beforeInterval;
				trajectory.entry = crossing.entry;
			} catch (const std::exception& e) {
				trajectory.status = cancel.stop_requested() ? NominalPropagationStatus::Cancelled
				                                            : NominalPropagationStatus::PropagationFailed;
				trajectory.error = e.what();
			}
			result.nominal = std::move(trajectory);
		}
		const auto earth = eph->queryMjd(EphemerisPlanet::Earth, result.epochMjdTdb);
		total += static_cast<int>(result.nodes.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(config.workers)
		for (int i = 0; i < static_cast<int>(result.nodes.size()); ++i) {
			auto& node = result.nodes[i];
			if (cancel.stop_requested()) {
				if (eligible(node))
					node.status = MovNodeStatus::Cancelled;
				continue;
			}
			if (eligible(node)) {
				node.weight = std::exp(-node.chi * node.chi / 2) * node.jacobianMov * node.jacobianSampling;
				try {
					const auto& x = node.state;
					const Vec3D r = Vec3D{x[0], x[1], x[2]} - earth.R, v = Vec3D{x[3], x[4], x[5]} - earth.V;
					if (v.dot(v) / 2 - Constants::GM_EARTH_ALONE / r.magnitude() < 0)
						node.status = MovNodeStatus::EarthBound;
					else {
						McClone clone;
						clone.index = i;
						clone.epochMjdTdb = result.epochMjdTdb;
						clone.state = x;
						auto p = propagator;
						mc::scan(clone, scan, p, *eph, cancel);
						node.entry = clone.entry;
						node.entryBeforeInterval = clone.entryBeforeInterval;
						node.status = cancel.stop_requested() ? MovNodeStatus::Cancelled : MovNodeStatus::Complete;
					}
				} catch (const std::exception& e) {
					node.status = MovNodeStatus::PropagationFailed;
					node.error = e.what();
				}
			}
			std::lock_guard lock(mutex);
			++done;
			if (progress && !callbackError)
				try {
					progress(done, total);
				} catch (...) {
					callbackError = std::current_exception();
				}
		}
		if (callbackError)
			std::rethrow_exception(callbackError);
		result.cancelled = cancel.stop_requested();
		double totalWeight = 0, hits = 0, missing = 0, outer = 0;
		int failedFits = 0, outerNodes = 0, earthBound = 0, missingMagnitude = 0;
		for (const auto& n : result.nodes) {
			if (n.status == MovNodeStatus::FitFailed)
				++failedFits;
			if (std::isfinite(n.chiSquare) && !n.absoluteMagnitude)
				++missingMagnitude;
			totalWeight += n.weight;
			if (n.status == MovNodeStatus::Complete && n.entry)
				hits += n.weight;
			if (n.status == MovNodeStatus::PropagationFailed || n.status == MovNodeStatus::Cancelled)
				missing += n.weight;
			if (n.status == MovNodeStatus::EarthBound)
				++earthBound;
			if (n.outerBoundary && n.weight > 0) {
				outer += n.weight;
				++outerNodes;
			}
		}
		if (!result.cancelled && totalWeight > 0) {
			result.probabilityLower = hits / totalWeight;
			result.probabilityUpper = (hits + missing) / totalWeight;
			result.outerBoundaryWeight = outer / totalWeight;
		}
		if (failedFits)
			result.warnings.push_back(std::to_string(failedFits) +
			                          " constrained fits failed; probability is conditional on resolved MOV nodes.");
		if (missingMagnitude)
			result.warnings.push_back(
			    std::to_string(missingMagnitude) +
			    " nodes lack usable photometry and were excluded by the absolute-magnitude selection.");
		if (outerNodes)
			result.warnings.push_back(
			    "Accepted nodes reach the sampled boundary; expand the domain and check grid convergence.");
		if (earthBound)
			result.warnings.push_back(std::to_string(earthBound) +
			                          " Earth-bound nodes retain normalization weight but are not propagated (NEO "
			                          "systematic-ranging policy).");
		if (missing > 0)
			result.warnings.push_back(
			    "Propagation failures have nonzero weight; only a probability bracket is available.");
		if (!totalWeight && !result.cancelled)
			throw std::runtime_error("MOV has no eligible probability mass");
		return result;
	}

}  // namespace orbital_core::hazard
