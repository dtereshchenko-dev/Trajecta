#include "orbital_core/optimizers/EigenLevenbergMarquardt.h"
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>
#include "logger/Logger.h"

#include <unsupported/Eigen/LevenbergMarquardt>

namespace orbital_core::optimizers {
	namespace {

		using Status = Eigen::LevenbergMarquardtSpace::Status;

		class ResidualFunctor final : public Eigen::DenseFunctor<double> {
		   public:
			ResidualFunctor(ILSQProblem& problem, std::function<bool()> cancel, const std::vector<double>& initial,
			                int maxEvaluations)
			    : Eigen::DenseFunctor<double>(problem.numParams(), 0),
			      m_problem(problem),
			      m_cancel(std::move(cancel)),
			      m_maxEvaluations(maxEvaluations)
			{
				if (cancelled() || !m_problem.acceptParams(initial))
					return;

				m_params = initial;
				m_problem.residuals(initial, m_residuals);
				++residualEvaluations;
				if (cancelled())
					return;

				if (m_residuals.size() != static_cast<std::size_t>(problem.numResiduals()))
					return;

				m_values = static_cast<int>(m_residuals.size());

				const Eigen::Map<const Eigen::VectorXd> r(m_residuals.data(), m_values);
				m_rejectionNorm = 2.0 * std::max(1.0, r.stableNorm());
				ready = m_values >= inputs() && r.allFinite() && std::isfinite(m_rejectionNorm);
			}

			int values() const { return m_values; }

			int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& residuals)
			{
				if (cancelled())
					return -1;

				lastTrialInvalid = false;
				const std::vector<double> params(x.data(), x.data() + x.size());
				if (!x.allFinite() || !m_problem.acceptParams(params))
					return reject(residuals);

				if (params != m_params) {
					if (residualEvaluations >= m_maxEvaluations) {
						evaluationLimited = true;
						return -1;
					}
					m_problem.residuals(params, m_residuals);
					++residualEvaluations;
					if (cancelled())
						return -1;

					if (m_residuals.size() != static_cast<std::size_t>(m_problem.numResiduals())) {
						failed = true;
						return -1;
					}

					m_params = params;
				}

				if (m_residuals.size() != static_cast<std::size_t>(m_values)) {
					failed = true;
					return -1;
				}

				residuals = Eigen::Map<const Eigen::VectorXd>(m_residuals.data(), m_values);
				if (!residuals.allFinite() || !std::isfinite(residuals.stableNorm()))
					return reject(residuals);

				return cancelled() ? -1 : 0;
			}

			int df(const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian)
			{
				Eigen::VectorXd residuals;
				if ((*this)(x, residuals) < 0 || lastTrialInvalid) {
					failed = !wasCancelled;
					return -1;
				}
				m_baseParams = m_params;
				m_rejectionsSinceBase = 0;
				m_problem.jacobian(m_params, jacobian);
				++jacobianEvaluations;
				if (cancelled())
					return -1;

				if (jacobian.rows() != m_values || jacobian.cols() != inputs() || !jacobian.allFinite()) {
					failed = true;
					return -1;
				}

				return cancelled() ? -1 : 0;
			}

			bool cancelled()
			{
				wasCancelled = wasCancelled || (m_cancel && m_cancel());
				return wasCancelled;
			}

			bool stoppedAtRejectedBase(const Eigen::VectorXd& x) const
			{
				return m_rejectionsSinceBase > 0 && std::equal(m_baseParams.begin(), m_baseParams.end(), x.data());
			}

			bool ready = false;
			bool failed = false;
			bool wasCancelled = false;
			bool lastTrialInvalid = false;
			bool evaluationLimited = false;
			int residualEvaluations = 0;
			int jacobianEvaluations = 0;
			int rejectedTrials = 0;

		   private:
			int reject(Eigen::VectorXd& residuals)
			{
				residuals = Eigen::VectorXd::Zero(m_values);
				residuals[0] = m_rejectionNorm;
				lastTrialInvalid = true;
				++rejectedTrials;
				++m_rejectionsSinceBase;
				return 0;
			}

			ILSQProblem& m_problem;
			std::function<bool()> m_cancel;
			std::vector<double> m_params;
			std::vector<double> m_residuals;
			std::vector<double> m_baseParams;
			int m_rejectionsSinceBase = 0;
			int m_maxEvaluations;
			int m_values = 0;
			double m_rejectionNorm = 0.0;
		};

		bool successful(Status status)
		{
			using namespace Eigen::LevenbergMarquardtSpace;
			return status == RelativeReductionTooSmall || status == RelativeErrorTooSmall ||
			       status == RelativeErrorAndReductionTooSmall || status == CosinusTooSmall;
		}

		const char* statusName(Status status)
		{
			using namespace Eigen::LevenbergMarquardtSpace;
			switch (status) {
				case NotStarted:
					return "not-started";
				case Running:
					return "iteration-limit";
				case ImproperInputParameters:
					return "invalid-input";
				case RelativeReductionTooSmall:
					return "function-tolerance";
				case RelativeErrorTooSmall:
					return "parameter-tolerance";
				case RelativeErrorAndReductionTooSmall:
					return "function-and-parameter-tolerance";
				case CosinusTooSmall:
					return "gradient-tolerance";
				case TooManyFunctionEvaluation:
					return "evaluation-limit";
				case FtolTooSmall:
					return "function-precision-limit";
				case XtolTooSmall:
					return "parameter-precision-limit";
				case GtolTooSmall:
					return "gradient-precision-limit";
				case UserAsked:
					return "callback-stop";
				default:
					return "unknown";
			}
		}

	}  // namespace

	constexpr const int maxIterations = 100;
	constexpr const int maxEvaluations = 1000;
	constexpr const double functionTolerance = 1e-12;
	constexpr const double parameterTolerance = 1e-11;
	constexpr const double gradientTolerance = 1e-12;
	constexpr const double initialStepFactor = 100.0;

	EigenLevenbergMarquardt::Result EigenLevenbergMarquardt::optimize(const std::vector<double>& initialParams,
	                                                                  ILSQProblem& problem,
	                                                                  std::function<bool()> isCancelled)
	{
		return optimize(initialParams, problem, std::move(isCancelled), {});
	}

	EigenLevenbergMarquardt::Result EigenLevenbergMarquardt::optimize(
	    const std::vector<double>& initialParams, ILSQProblem& problem, std::function<bool()> isCancelled,
	    std::function<bool(const std::vector<double>&)> updateSelection)
	{
		if (initialParams.empty() || initialParams.size() != static_cast<std::size_t>(problem.numParams()) ||
		    problem.numResiduals() <= 0 ||
		    !std::all_of(initialParams.begin(), initialParams.end(), [](double v) { return std::isfinite(v); }))
			throw std::invalid_argument("Eigen LM: invalid parameters or configuration");

		Result result;
		result.params = initialParams;
		result.finalCost = std::numeric_limits<double>::quiet_NaN();

		Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(initialParams.data(), initialParams.size());
		bool restart;
		do {
			restart = false;
			std::vector<double> lastChecked;
			result.converged = false;
			result.finalCost = std::numeric_limits<double>::quiet_NaN();
			ResidualFunctor functor(problem, isCancelled, result.params, maxEvaluations);
			if (!functor.ready) {
				result.cancelled = functor.wasCancelled;
				return result;
			}
			Eigen::LevenbergMarquardt<ResidualFunctor> solver(functor);
			solver.setFtol(functionTolerance);
			solver.setXtol(parameterTolerance);
			solver.setGtol(gradientTolerance);
			solver.setFactor(initialStepFactor);
			solver.setMaxfev(maxEvaluations);
			Status status = solver.minimizeInit(x);
			LOG_INFO << "[EigenLM] start n=" << x.size() << " m=" << functor.values() << " ftol=" << functionTolerance
			         << " xtol=" << parameterTolerance << " gtol=" << gradientTolerance
			         << " factor=" << initialStepFactor;
			for (int iteration = 0;
			     iteration < maxIterations && (status == Eigen::LevenbergMarquardtSpace::NotStarted ||
			                                   status == Eigen::LevenbergMarquardtSpace::Running);
			     ++iteration) {
				if (functor.cancelled())
					break;
				const Eigen::VectorXd previous = x;
				status = solver.minimizeOneStep(x);
				++result.iterations;
				result.params.assign(x.data(), x.data() + x.size());
				if (functor.wasCancelled)
					break;
				if (updateSelection && !functor.failed && !functor.lastTrialInvalid && !functor.evaluationLimited &&
				    (status == Eigen::LevenbergMarquardtSpace::Running || successful(status)) &&
				    (x.array() != previous.array()).any()) {
					lastChecked = result.params;
					if (updateSelection(result.params)) {
						restart = true;
						break;
					}
				}
			}
			if (restart)
				continue;
			const bool validStop = !functor.cancelled() && !functor.failed && !functor.lastTrialInvalid &&
			                       !functor.evaluationLimited && !functor.stoppedAtRejectedBase(x);
			if (!functor.cancelled())
				result.finalCost = problem.displayCost(result.params);
			result.cancelled = functor.cancelled();
			result.converged = validStop && !result.cancelled && std::isfinite(result.finalCost) && successful(status);
			if (result.converged && updateSelection && lastChecked != result.params && updateSelection(result.params)) {
				restart = true;
				continue;
			}
			LOG_INFO << "[EigenLM] finish status="
			         << (functor.evaluationLimited ? "evaluation-limit" : statusName(status))
			         << " converged=" << result.converged << " cancelled=" << result.cancelled
			         << " iterations=" << result.iterations << " f=" << functor.residualEvaluations
			         << " J=" << functor.jacobianEvaluations << " invalidTrials=" << functor.rejectedTrials;
		} while (restart);

		return result;
	}

}  // namespace orbital_core::optimizers
