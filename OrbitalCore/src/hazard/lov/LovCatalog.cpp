#include "hazard/lov/LovTraceDetail.h"
#include <utility>
#include <vector>

#include "hazard/ObservationNormalMetric.h"
#include "hazard/HazardAudit.h"
#include "logger/Logger.h"

#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace orbital_core::hazard::lov_detail {
	namespace {

		constexpr int kCorrectorIterations = 20;
		constexpr double kCorrectionTolerance = 1.0e-3;

		struct ScaledWeakField {
			VecX units;
			VecX unitDirection;
			VecX coordinateStepPerSigma;
			MatX coordinateCovariance;
			MatX physicalCovariance;
			VecX residuals;
			MatX jacobianScaled;
			double weakSigma{};
			double chiSquare{};
			bool valid{false};
		};

		struct CatalogFailureInfo {
			const char* stage{"unknown"};
			int iteration{-1};
			double correctionNorm{};
			double rms{};
		};

		void auditField(const char* stage, const VecX& coordinates, const ScaledWeakField& field, int iteration = -1)
		{
			audit::write(stage, [&](auto& out) {
				out << " iteration=" << iteration << " valid=" << field.valid << " weakSigma=" << field.weakSigma
				    << " chiSquare=" << field.chiSquare;
				for (Eigen::Index i = 0; i < coordinates.size(); ++i) {
					out << " equ" << i << '=' << coordinates[i];
					if (!field.valid)
						continue;
					out << " unit" << i << '=' << field.units[i] << " dir" << i << '=' << field.unitDirection[i];
					for (Eigen::Index j = 0; j < coordinates.size(); ++j)
						out << " cov" << i << j << '=' << field.coordinateCovariance(i, j);
				}
			});
		}

		void alignDirection(VecX& direction, const VecX* referenceDirection)
		{
			if (referenceDirection && referenceDirection->size() == direction.size() &&
			    direction.dot(*referenceDirection) < 0.0) {
				direction = -direction;
				return;
			}

			if (!referenceDirection && direction.size() > 0 && direction(0) < 0.0) {
				direction = -direction;
			}
		}

		[[nodiscard]] ScaledWeakField fieldFromCovariance(const TraceEvaluationContext& ctx, const VecX& coordinates,
		                                                  const FitParameterCovariance& coordinateCovariance,
		                                                  const FitParameterCovariance& physicalCovariance,
		                                                  const VecX* referenceDirection = nullptr)
		{
			ScaledWeakField out;
			const VecX units = lovCoordinateUnits(ctx.chart.coordinateSystem, coordinates, ctx.cfg.coordinateScaling);
			if (units.size() != coordinates.size()) {
				return out;
			}

			const VecX scales = units.cwiseInverse();
			const MatX scaledCovariance = scales.asDiagonal() * coordinateCovariance * scales.asDiagonal();
			LovWeakDirection weak = weakDirectionFromCovariance(scaledCovariance);
			if (!weak.valid || weak.parameterStepPerSigma.size() != coordinates.size() ||
			    !(weak.covarianceEigenvalue > 0.0)) {
				return out;
			}

			out.unitDirection = weak.parameterStepPerSigma / std::sqrt(weak.covarianceEigenvalue);
			alignDirection(out.unitDirection, referenceDirection);
			out.units = units;
			out.coordinateStepPerSigma = units.asDiagonal() * out.unitDirection * std::sqrt(weak.covarianceEigenvalue);
			out.weakSigma = std::sqrt(weak.covarianceEigenvalue);
			out.coordinateCovariance = coordinateCovariance;
			out.physicalCovariance = physicalCovariance;
			out.valid = out.coordinateStepPerSigma.array().isFinite().all();
			return out;
		}

		[[nodiscard]] ScaledWeakField fieldFromObservations(const TraceEvaluationContext& ctx, const VecX& coordinates,
		                                                    const VecX* referenceDirection = nullptr)
		{
			ScaledWeakField out;
			const auto& linearization = ctx.cfg.impactRefinement.rmsScaling.linearization;
			if (!linearization)
				return out;

			const auto physical = physicalParametersFromChart(ctx.chart, coordinates);
			const auto physicalJacobian = physicalJacobianFromChart(ctx.chart, coordinates);
			if (!physical || !physicalJacobian)
				return out;

			const auto lin = linearization(*physical);
			const Eigen::Index n = coordinates.size();
			if (!lin.valid || lin.jacobian.cols() != n ||
			    lin.jacobian.rows() != static_cast<Eigen::Index>(lin.residuals.size()) || lin.residuals.empty() ||
			    !lin.jacobian.array().isFinite().all()) {
				return out;
			}

			const VecX units = lovCoordinateUnits(ctx.chart.coordinateSystem, coordinates, ctx.cfg.coordinateScaling);
			if (units.size() != coordinates.size())
				return out;

			out.residuals =
			    Eigen::Map<const VecX>(lin.residuals.data(), static_cast<Eigen::Index>(lin.residuals.size()));
			out.jacobianScaled = lin.jacobian * (*physicalJacobian) * units.asDiagonal();

			const ObservationNormalMetric metric = makeObservationNormalMetric(out.jacobianScaled, lin.sumSquares);
			if (!metric.valid)
				return out;

			const VecX& values = metric.eigenvalues;
			if (values.size() != n || !(values(values.size() - 1) > 0.0))
				return out;

			Eigen::Index weakIndex = -1;
			for (Eigen::Index i = 0; i < values.size(); ++i) {
				if (std::isfinite(values(i)) && values(i) > metric.eigenvalueFloor) {
					weakIndex = i;
					break;
				}
			}
			if (weakIndex < 0)
				return out;

			out.unitDirection = metric.eigenvectors.col(weakIndex);
			alignDirection(out.unitDirection, referenceDirection);
			out.units = units;

			const double chiSquare = lin.sumSquares;
			const double scale = metric.covarianceScale;
			const double weakSigma = std::sqrt(scale / values(weakIndex));
			out.weakSigma = weakSigma;
			out.coordinateStepPerSigma = units.asDiagonal() * out.unitDirection * weakSigma;

			out.coordinateCovariance = units.asDiagonal() * metric.covariance * units.asDiagonal();
			out.physicalCovariance = (*physicalJacobian) * out.coordinateCovariance * physicalJacobian->transpose();
			out.physicalCovariance = 0.5 * (out.physicalCovariance + out.physicalCovariance.transpose()).eval();
			out.chiSquare = chiSquare;
			out.valid = out.residuals.array().isFinite().all() && out.jacobianScaled.array().isFinite().all() &&
			            out.coordinateStepPerSigma.array().isFinite().all() &&
			            out.physicalCovariance.array().isFinite().all();
			return out;
		}

		[[nodiscard]] MatX orthogonalComplement(const VecX& unitDirection)
		{
			const Eigen::Index n = unitDirection.size();
			MatX basis = MatX::Zero(n, std::max<Eigen::Index>(0, n - 1));
			Eigen::Index count = 0;

			for (Eigen::Index axis = 0; axis < n && count < n - 1; ++axis) {
				VecX v = VecX::Unit(n, axis);
				v -= unitDirection * unitDirection.dot(v);
				for (Eigen::Index j = 0; j < count; ++j)
					v -= basis.col(j) * basis.col(j).dot(v);

				const double norm = v.norm();
				if (norm <= 1.0e-10)
					continue;
				basis.col(count++) = v / norm;
			}

			if (count != n - 1)
				return {};
			return basis;
		}

		[[nodiscard]] LovCatalogPoint pointFromField(const TraceEvaluationContext& ctx, int sampleIndex, double sigma,
		                                             const VecX& coordinates, const ScaledWeakField& field)
		{
			LovCatalogPoint point;
			point.sampleIndex = sampleIndex;
			point.sigma = sigma;
			point.coordinates = coordinates;
			point.coordinateCovariance = field.coordinateCovariance;
			point.physicalCovariance = field.physicalCovariance;
			point.chiSquare = field.chiSquare;

			const auto physical = physicalParametersFromChart(ctx.chart, coordinates);
			if (!physical)
				return point;
			point.parameters = *physical;

			const auto physicalStep =
			    physicalStepPerSigmaFromChart(ctx.chart, coordinates, field.coordinateStepPerSigma);
			if (!physicalStep)
				return point;
			point.physicalStepPerSigma = *physicalStep;
			point.valid = field.valid && point.parameters.array().isFinite().all() &&
			              point.physicalStepPerSigma.array().isFinite().all();
			return point;
		}

		[[nodiscard]] LovCatalogPoint linearCatalogPoint(const TraceEvaluationContext& ctx, int sampleIndex,
		                                                 double sigma, double coordinateSigma)
		{
			const ScaledWeakField weak =
			    fieldFromCovariance(ctx, ctx.chart.nominalCoordinates, ctx.chart.coordinateCovariance, ctx.covariance);
			if (!weak.valid)
				return {};

			const VecX coordinates = ctx.chart.nominalCoordinates + coordinateSigma * weak.coordinateStepPerSigma;
			return pointFromField(ctx, sampleIndex, sigma, coordinates, weak);
		}

		[[nodiscard]] double linearCoordinateSigma(const std::vector<double>& sigmas, std::size_t index,
		                                           std::size_t center) noexcept
		{
			if (index == center || center + 1 >= sigmas.size())
				return sigmas[index];

			const double firstDelta = sigmas[center + 1] - sigmas[center];
			if (index < center) {
				const double outerDelta =
				    index > 0 ? sigmas[index] - sigmas[index - 1] : sigmas[index + 1] - sigmas[index];
				return sigmas[index] - (outerDelta - firstDelta);
			}

			const double outerDelta =
			    index + 1 < sigmas.size() ? sigmas[index + 1] - sigmas[index] : sigmas[index] - sigmas[index - 1];
			return sigmas[index] + (outerDelta - firstDelta);
		}

		[[nodiscard]] LovCatalogPoint constrainedCatalogPoint(const TraceEvaluationContext& ctx, int sampleIndex,
		                                                      double sigma, VecX coordinates,
		                                                      VecX& continuationDirection,
		                                                      CatalogFailureInfo* failure = nullptr)
		{
			VecX direction = continuationDirection;
			double correctionNorm = 0.0;
			bool converged = false;
			double previousRms = std::numeric_limits<double>::infinity();
			int stalledIterations = 0;
			bool stalledConvergence = false;
			ScaledWeakField convergedField;

			for (int iteration = 0; iteration < kCorrectorIterations; ++iteration) {
				ScaledWeakField field = fieldFromObservations(ctx, coordinates, &direction);
				auditField("catalog-corrector-field", coordinates, field, iteration);
				if (!field.valid) {
					if (failure) {
						failure->stage = "corrector-field";
						failure->iteration = iteration;
					}
					return {};
				}
				direction = field.unitDirection;


				const MatX complement = orthogonalComplement(direction);
				if (complement.cols() + 1 != coordinates.size()) {
					if (failure) {
						failure->stage = "corrector-basis";
						failure->iteration = iteration;
					}
					return {};
				}

				const MatX reducedJacobian = field.jacobianScaled * complement;
				const VecX reducedCorrection =
				    reducedJacobian.completeOrthogonalDecomposition().solve(-field.residuals);
				if (reducedCorrection.size() != complement.cols() || !reducedCorrection.array().isFinite().all()) {
					if (failure) {
						failure->stage = "corrector-solve";
						failure->iteration = iteration;
					}
					return {};
				}

				const VecX scaledCorrection = complement * reducedCorrection;
				correctionNorm = (field.jacobianScaled * scaledCorrection).norm() /
				                 std::sqrt(static_cast<double>(coordinates.size()));
				const double rms =
				    std::sqrt(field.residuals.squaredNorm() / static_cast<double>(field.residuals.size()));
				if (failure) {
					failure->iteration = iteration;
					failure->correctionNorm = correctionNorm;
					failure->rms = rms;
				}

				VecX coordinateCorrection = field.units.asDiagonal() * scaledCorrection;
				audit::write("catalog-correction", [&](auto& out) {
					out << " iteration=" << iteration << " correctionNorm=" << correctionNorm << " rms=" << rms;
					for (Eigen::Index i = 0; i < coordinates.size(); ++i)
						out << " delta" << i << '=' << coordinateCorrection[i];
				});
				bool accepted = false;
				for (int relaxation = 0; relaxation < 4; ++relaxation) {
					const VecX candidate = coordinates + coordinateCorrection;
					if (physicalParametersFromChart(ctx.chart, candidate)) {
						coordinates = candidate;
						accepted = true;
						break;
					}
					coordinateCorrection *= 0.5;
				}
				if (!accepted) {
					const VecX candidate = coordinates + coordinateCorrection;
					if (physicalParametersFromChart(ctx.chart, candidate)) {
						coordinates = candidate;
					} else {
						if (failure)
							failure->stage = "corrector-orbit";
						return {};
					}
				}

				if (correctionNorm < kCorrectionTolerance) {
					converged = true;
					convergedField = std::move(field);
					break;
				}

				if (std::isfinite(previousRms)) {
					if (rms > previousRms * 1.1) {
						++stalledIterations;
						stalledConvergence = false;
					} else if (rms > previousRms * 0.999) {
						++stalledIterations;
						stalledConvergence = true;
					}
					if (stalledIterations > 5) {
						if (!stalledConvergence) {
							if (failure)
								failure->stage = "corrector-increasing";
							return {};
						}
						converged = true;
						convergedField = std::move(field);
						break;
					}
				}
				previousRms = rms;
			}

			if (!converged) {
				if (failure)
					failure->stage = "corrector-nonconvergent";
				return {};
			}

			continuationDirection = convergedField.unitDirection;

			const auto physicalJacobian = physicalJacobianFromChart(ctx.chart, coordinates);
			if (!physicalJacobian) {
				if (failure)
					failure->stage = "corrector-final-jacobian";
				return {};
			}

			FitParameterCovariance physicalCovariance =
			    ((*physicalJacobian) * convergedField.coordinateCovariance * physicalJacobian->transpose()).eval();
			physicalCovariance = 0.5 * (physicalCovariance + physicalCovariance.transpose()).eval();
			ScaledWeakField storedField = fieldFromCovariance(ctx, coordinates, convergedField.coordinateCovariance,
			                                                  physicalCovariance, &continuationDirection);
			if (!storedField.valid) {
				if (failure)
					failure->stage = "corrector-final-field";
				return {};
			}
			storedField.chiSquare = convergedField.chiSquare;
			return pointFromField(ctx, sampleIndex, sigma, coordinates, storedField);
		}

		struct CatalogContinuationState {
			VecX direction;
			VecX units;
			double weakSigma{};
		};

		[[nodiscard]] std::optional<VecX> sigmaPropagationStep(const TraceEvaluationContext& ctx,
		                                                       const VecX& coordinates, double sigmaStep,
		                                                       CatalogContinuationState& state)
		{
			if (state.direction.size() != coordinates.size() || state.units.size() != coordinates.size() ||
			    !(state.weakSigma > 0.0)) {
				return std::nullopt;
			}

			const VecX midpoint =
			    coordinates + 0.5 * sigmaStep * state.units.asDiagonal() * state.direction * state.weakSigma;
			audit::write("catalog-predictor-input", [&](auto& out) {
				out << " step=" << sigmaStep << " weakSigma=" << state.weakSigma;
				for (Eigen::Index i = 0; i < coordinates.size(); ++i)
					out << " equ" << i << '=' << coordinates[i] << " unit" << i << '=' << state.units[i] << " dir" << i
					    << '=' << state.direction[i];
			});
			if (!physicalParametersFromChart(ctx.chart, midpoint))
				return std::nullopt;

			const ScaledWeakField midpointField = fieldFromObservations(ctx, midpoint, &state.direction);
			auditField("catalog-midpoint", midpoint, midpointField);
			if (!midpointField.valid)
				return std::nullopt;


			state.units = midpointField.units;
			const VecX predicted = coordinates + sigmaStep * state.units.asDiagonal() * midpointField.unitDirection *
			                                         midpointField.weakSigma;
			if (!physicalParametersFromChart(ctx.chart, predicted))
				return std::nullopt;
			return predicted;
		}

		[[nodiscard]] std::optional<VecX> propagateCoordinatesToSigma(const TraceEvaluationContext& ctx,
		                                                              const VecX& coordinates, double requestedStep,
		                                                              CatalogContinuationState& state,
		                                                              CatalogFailureInfo* failure)
		{
			const CatalogContinuationState entryState = state;
			double step = requestedStep;
			double completed = 0.0;
			VecX current = coordinates;

			for (;;) {
				const auto predicted = sigmaPropagationStep(ctx, current, step, state);
				if (!predicted) {
					if (std::abs(step) >= std::abs(requestedStep) / 130.0) {
						step *= 0.5;
						continue;
					}
					if (failure)
						failure->stage = "predictor-halving";
					return std::nullopt;
				}

				current = *predicted;
				completed += step;
				if (std::abs(completed) >= std::abs(requestedStep))
					return current;

				const ScaledWeakField partialField = fieldFromObservations(ctx, current, &entryState.direction);
				if (!partialField.valid) {
					if (failure)
						failure->stage = "predictor-partial-field";
					return std::nullopt;
				}
				state.direction = partialField.unitDirection;
				state.units = partialField.units;
				state.weakSigma = partialField.weakSigma;
			}
		}

		[[nodiscard]] LovCatalogPoint continueCatalogPoint(const TraceEvaluationContext& ctx,
		                                                   const LovCatalogPoint& previous, int sampleIndex,
		                                                   double sigma, double predictorStep,
		                                                   CatalogContinuationState& state,
		                                                   CatalogFailureInfo* failure = nullptr)
		{
			const audit::Scope auditScope(sampleIndex, sigma);
			const auto predicted =
			    propagateCoordinatesToSigma(ctx, previous.coordinates, predictorStep, state, failure);
			if (!predicted)
				return {};
			return constrainedCatalogPoint(ctx, sampleIndex, sigma, *predicted, state.direction, failure);
		}

	}  // namespace

	std::vector<LovCatalogPoint> buildLovCatalog(const TraceEvaluationContext& ctx, const std::vector<double>& sigmas)
	{
		std::vector<LovCatalogPoint> catalog(sigmas.size());
		const auto centerIt = std::find(sigmas.begin(), sigmas.end(), 0.0);
		if (centerIt == sigmas.end())
			return catalog;

		const std::size_t center = static_cast<std::size_t>(std::distance(sigmas.begin(), centerIt));

		if (!ctx.cfg.impactRefinement.rmsScaling.linearization) {
			for (std::size_t i = 0; i < sigmas.size(); ++i) {
				catalog[i] =
				    linearCatalogPoint(ctx, static_cast<int>(i), sigmas[i], linearCoordinateSigma(sigmas, i, center));
			}
			LOG_INFO << "[LineOfVariation] nonlinear catalog unavailable"
			         << "  fallback=linear"
			         << "  reason=no-observation-linearization";
			return catalog;
		}

		catalog[center] = linearCatalogPoint(ctx, static_cast<int>(center), 0.0, 0.0);
		if (!catalog[center].valid)
			return catalog;

		const double referenceChiSquare = nominalObservationChiSquareRaw(ctx);
		if (!std::isfinite(referenceChiSquare))
			return catalog;
		catalog[center].chiSquare = referenceChiSquare;

		const ScaledWeakField nominalField =
		    fieldFromCovariance(ctx, ctx.chart.nominalCoordinates, ctx.chart.coordinateCovariance, ctx.covariance);
		if (!nominalField.valid)
			return catalog;
		const CatalogContinuationState nominalState{nominalField.unitDirection, nominalField.units,
		                                            nominalField.weakSigma};

		CatalogContinuationState positiveState = nominalState;
		for (std::size_t i = center + 1; i < sigmas.size(); ++i) {
			CatalogFailureInfo failure;
			catalog[i] = continueCatalogPoint(ctx, catalog[i - 1], static_cast<int>(i), sigmas[i],
			                                  catalogPredictorStep(sigmas, i, center), positiveState, &failure);
			if (!catalog[i].valid) {
				LOG_INFO << "[LineOfVariation] nonlinear catalog stopped"
				         << "  branch=positive  sample=" << i << "  sigma=" << sigmas[i] << "  stage=" << failure.stage
				         << "  iteration=" << failure.iteration << "  correctionNorm=" << failure.correctionNorm
				         << "  rms=" << failure.rms;
				break;
			}

			const double sigmaQ = catalogSigmaQ(catalog[i].chiSquare, referenceChiSquare);
			if (std::isfinite(sigmaQ) && sigmaQ > ctx.cfg.sigmaMax) {
				LOG_INFO << "[LineOfVariation] nonlinear catalog stopped"
				         << "  branch=positive  sample=" << i << "  sigma=" << sigmas[i] << "  stage=sigma-q-limit"
				         << "  sigmaQ=" << sigmaQ;
				break;
			}
		}

		CatalogContinuationState negativeState = nominalState;
		for (std::size_t i = center; i-- > 0;) {
			CatalogFailureInfo failure;
			catalog[i] = continueCatalogPoint(ctx, catalog[i + 1], static_cast<int>(i), sigmas[i],
			                                  catalogPredictorStep(sigmas, i, center), negativeState, &failure);
			if (!catalog[i].valid) {
				LOG_INFO << "[LineOfVariation] nonlinear catalog stopped"
				         << "  branch=negative  sample=" << i << "  sigma=" << sigmas[i] << "  stage=" << failure.stage
				         << "  iteration=" << failure.iteration << "  correctionNorm=" << failure.correctionNorm
				         << "  rms=" << failure.rms;
				break;
			}

			const double sigmaQ = catalogSigmaQ(catalog[i].chiSquare, referenceChiSquare);
			if (std::isfinite(sigmaQ) && sigmaQ > ctx.cfg.sigmaMax) {
				LOG_INFO << "[LineOfVariation] nonlinear catalog stopped"
				         << "  branch=negative  sample=" << i << "  sigma=" << sigmas[i] << "  stage=sigma-q-limit"
				         << "  sigmaQ=" << sigmaQ;
				break;
			}
		}
		return catalog;
	}

	std::optional<LovCatalogPoint> catalogPointAtSigma(const TraceEvaluationContext& ctx, double sigma)
	{
		if (!ctx.catalog || ctx.catalog->empty())
			return linearCatalogPoint(ctx, -1, sigma, sigma);

		const LovCatalogPoint* nearest = catalogAnchorAtSigma(*ctx.catalog, sigma);
		if (!nearest)
			return std::nullopt;
		const double nearestDistance = std::abs(nearest->sigma - sigma);
		if (nearestDistance <= 1.0e-14)
			return *nearest;

		if (!ctx.cfg.impactRefinement.rmsScaling.linearization) {
			const LovCatalogPoint* lower = nullptr;
			const LovCatalogPoint* upper = nullptr;
			for (const LovCatalogPoint& point : *ctx.catalog) {
				if (!point.valid)
					continue;
				if (point.sigma < sigma && (!lower || point.sigma > lower->sigma)) {
					lower = &point;
				}
				if (point.sigma > sigma && (!upper || point.sigma < upper->sigma)) {
					upper = &point;
				}
			}
			if (!lower || !upper || !(upper->sigma > lower->sigma))
				return *nearest;

			const double fraction = (sigma - lower->sigma) / (upper->sigma - lower->sigma);
			const VecX coordinates = lower->coordinates + fraction * (upper->coordinates - lower->coordinates);
			ScaledWeakField field =
			    fieldFromCovariance(ctx, ctx.chart.nominalCoordinates, ctx.chart.coordinateCovariance, ctx.covariance);
			if (!field.valid)
				return std::nullopt;
			field.coordinateStepPerSigma = (upper->coordinates - lower->coordinates) / (upper->sigma - lower->sigma);
			return pointFromField(ctx, -1, sigma, coordinates, field);
		}

		const LovCatalogPoint* previous = nullptr;
		const LovCatalogPoint* next = nullptr;
		for (const LovCatalogPoint& candidate : *ctx.catalog) {
			if (!candidate.valid)
				continue;
			if (candidate.sigma < nearest->sigma && (!previous || candidate.sigma > previous->sigma)) {
				previous = &candidate;
			}
			if (candidate.sigma > nearest->sigma && (!next || candidate.sigma < next->sigma)) {
				next = &candidate;
			}
		}

		ScaledWeakField nearestField =
		    fieldFromCovariance(ctx, nearest->coordinates, nearest->coordinateCovariance, nearest->physicalCovariance);
		if (!nearestField.valid)
			return std::nullopt;

		const VecX* neighborDifference = nullptr;
		VecX difference;
		if (next) {
			difference = next->coordinates - nearest->coordinates;
			neighborDifference = &difference;
		} else if (previous) {
			difference = nearest->coordinates - previous->coordinates;
			neighborDifference = &difference;
		}
		if (neighborDifference && neighborDifference->dot(nearestField.unitDirection) < 0.0) {
			nearestField.unitDirection = -nearestField.unitDirection;
		}

		CatalogContinuationState interpolationState{nearestField.unitDirection, nearestField.units,
		                                            nearestField.weakSigma};
		LovCatalogPoint point =
		    continueCatalogPoint(ctx, *nearest, -1, sigma, sigma - nearest->sigma, interpolationState);
		if (!point.valid)
			return std::nullopt;
		return point;
	}

}  // namespace orbital_core::hazard::lov_detail
