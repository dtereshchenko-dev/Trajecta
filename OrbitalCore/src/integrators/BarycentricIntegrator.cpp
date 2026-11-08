#include "BarycentricIntegrator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include <limits>
#include <utility>

namespace {
	class BarycentricDynamics final : public IIntegrationDynamics {
	   public:
		BarycentricDynamics(const IIntegrationDynamics& dynamics, const DE440Source& source)
		    : m_dynamics(dynamics), m_source(source)
		{
		}

		const EphemerisForceTerms& sun(double epoch, double offset = 0.0) const
		{
			if (epoch != m_epoch || offset != m_offset) {
				Vec3D positions[8], velocities[8];
				m_source.getAllPlanetStatesOffsetMjd(epoch, offset, positions, velocities, nullptr, nullptr, &m_sun);
				m_epoch = epoch;
				m_offset = offset;
			}
			return m_sun;
		}

		template <class State>
		State heliocentric(State state, double epoch, double offset = 0.0) const
		{
			const auto& terms = sun(epoch, offset);
			state.R -= terms.sunPosition;
			state.V -= terms.sunVelocity;
			return state;
		}

		template <class State>
		State barycentric(State state, double mjd) const
		{
			const auto& terms = sun(mjd);
			state.R += terms.sunPosition;
			state.V += terms.sunVelocity;
			return state;
		}

		Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng) const override
		{
			return accelerationAtOffsetMjd(r, v, mjd, 0.0, ng);
		}

		Vec3D accelerationAtOffsetMjd(const Vec3D& r, const Vec3D& v, double epoch, double offset,
		                              const NonGravModel& ng) const override
		{
			const auto terms = sun(epoch, offset);
			return m_dynamics.accelerationAtOffsetMjd(r - terms.sunPosition, v - terms.sunVelocity, epoch, offset, ng) +
			       terms.sunAcceleration;
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng, Vec3D& a,
		                                Eigen::Matrix3d& ar, Eigen::Matrix3d& av) const override
		{
			accelerationAndPartialsOffsetMjd(r, v, mjd, 0.0, ng, a, ar, av);
		}

		void accelerationAndPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epoch, double offset,
		                                      const NonGravModel& ng, Vec3D& a, Eigen::Matrix3d& ar,
		                                      Eigen::Matrix3d& av) const override
		{
			const auto terms = sun(epoch, offset);
			m_dynamics.accelerationAndPartialsOffsetMjd(r - terms.sunPosition, v - terms.sunVelocity, epoch, offset, ng,
			                                            a, ar, av);
			a += terms.sunAcceleration;
		}

		Eigen::Matrix<double, 3, 4> nonGravPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd,
		                                               const NonGravModel& ng) const override
		{
			return nonGravPartialsOffsetMjd(r, v, mjd, 0.0, ng);
		}

		Eigen::Matrix<double, 3, 4> nonGravPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epoch,
		                                                     double offset, const NonGravModel& ng) const override
		{
			const auto terms = sun(epoch, offset);
			return m_dynamics.nonGravPartialsOffsetMjd(r - terms.sunPosition, v - terms.sunVelocity, epoch, offset, ng);
		}

		IntegrationOptions options(IntegrationOptions value) const
		{
			if (value.nextEncounterStep)
				value.nextEncounterStep = [this, next = value.nextEncounterStep](double step, double mjd,
				                                                                 const IntegratedState& state) mutable {
					return next(step, mjd, heliocentric(state, mjd));
				};
			if (value.nextEncounterStepOffset)
				value.nextEncounterStepOffset = [this, next = value.nextEncounterStepOffset](
				                                    double step, double epoch, double offset,
				                                    const IntegratedState& state) mutable {
					return next(step, epoch, offset, heliocentric(state, epoch, offset));
				};
			if (value.variationalObserver)
				value.variationalObserver = [this, observe = value.variationalObserver](
				                                double mjd, const IntegratedStateFull& state) mutable {
					return observe(mjd, heliocentric(state, mjd));
				};
			if (value.variationalObserverOffset)
				value.variationalObserverOffset = [this, observe = value.variationalObserverOffset](
				                                      double epoch, double offset,
				                                      const IntegratedStateFull& state) mutable {
					return observe(epoch, offset, heliocentric(state, epoch, offset));
				};
			if (value.denseStepObserver)
				value.denseStepObserver = [source = &m_source, observe = value.denseStepObserver](DenseIntegratedStep step) {
					auto evaluate = std::move(step.evaluate);
					step.evaluate = [source, evaluate = std::move(evaluate)](double mjd) {
						auto state = evaluate(mjd);
						Vec3D positions[8], velocities[8];
						EphemerisForceTerms terms;
						source->getAllPlanetStatesMjd(mjd, positions, velocities, nullptr, nullptr, &terms);
						state.R -= terms.sunPosition;
						state.V -= terms.sunVelocity;
						return state;
					};
					observe(std::move(step));
				};
			return value;
		}

		template <class State>
		IntegrationReport finish(IntegrationReport report, const State& initial, double epoch,
		                         std::span<const double> targets, std::span<State> outputs) const
		{
			for (std::size_t i = 0; i < report.completedOutputs; ++i)
				outputs[i] = targets[i] == epoch ? initial : heliocentric(outputs[i], targets[i]);
			report.finalState = report.reachedMjd == epoch ? IntegratedState{initial.R, initial.V}
			                                               : heliocentric(report.finalState, report.reachedMjd);
			return report;
		}

	   private:
		const IIntegrationDynamics& m_dynamics;
		const DE440Source& m_source;
		mutable double m_epoch = std::numeric_limits<double>::quiet_NaN();
		mutable double m_offset = 0.0;
		mutable EphemerisForceTerms m_sun;
	};

	const DE440Source* availableSource()
	{
		const auto* source = dynamic_cast<DE440Source*>(getDefaultEphemerisSource());
		return source && source->isAvailable() ? source : nullptr;
	}
}  // namespace

BarycentricIntegrator::BarycentricIntegrator(IIntegrator_p engine) : m_engine(std::move(engine)) {}

const char* BarycentricIntegrator::name() const noexcept { return m_engine->name(); }

IntegrationReport BarycentricIntegrator::integrateState(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
                                                        const IntegratedState& initial, double epoch,
                                                        std::span<const double> targets,
                                                        std::span<IntegratedState> outputs,
                                                        const IntegrationOptions& options,
                                                        const AcceptedStateStepObserver& observer) const
{
	const auto* source = availableSource();
	if (!source || options.cancel.stop_requested() || targets.empty())
		return m_engine->integrateState(dynamics, ng, initial, epoch, targets, outputs, options, observer);
	const BarycentricDynamics translated(dynamics, *source);
	AcceptedStateStepObserver translatedObserver;
	if (observer)
		translatedObserver = [&translated, &observer](const AcceptedStateStep& step) {
			AcceptedStateStep translatedStep{step.beginMjd, step.endMjd,
			                                 translated.heliocentric(step.begin, step.beginMjd),
			                                 translated.heliocentric(step.end, step.endMjd)};
			if (step.interpolateMjd)
				translatedStep.interpolateMjd = [&translated, &step](double mjd) {
					return translated.heliocentric(step.interpolateMjd(mjd), mjd);
				};
			return observer(translatedStep);
		};
	const auto report = m_engine->integrateState(translated, ng, translated.barycentric(initial, epoch), epoch, targets,
	                                             outputs, translated.options(options), translatedObserver);
	return translated.finish(report, initial, epoch, targets, outputs);
}

IntegrationReport BarycentricIntegrator::integrateStateWithStm(
    const IIntegrationDynamics& dynamics, const NonGravModel& ng, const IntegratedStateWithStm& initial, double epoch,
    std::span<const double> targets, std::span<IntegratedStateWithStm> outputs, const IntegrationOptions& options) const
{
	const auto* source = availableSource();
	if (!source || options.cancel.stop_requested() || targets.empty())
		return m_engine->integrateStateWithStm(dynamics, ng, initial, epoch, targets, outputs, options);
	const BarycentricDynamics translated(dynamics, *source);
	const auto report = m_engine->integrateStateWithStm(translated, ng, translated.barycentric(initial, epoch), epoch,
	                                                    targets, outputs, translated.options(options));
	return translated.finish(report, initial, epoch, targets, outputs);
}

IntegrationReport BarycentricIntegrator::integrateStateWithStmAndSens(
    const IIntegrationDynamics& dynamics, const NonGravModel& ng, const IntegratedStateFull& initial, double epoch,
    std::span<const double> targets, std::span<IntegratedStateFull> outputs, const IntegrationOptions& options) const
{
	const auto* source = availableSource();
	if (!source || options.cancel.stop_requested() || targets.empty())
		return m_engine->integrateStateWithStmAndSens(dynamics, ng, initial, epoch, targets, outputs, options);
	const BarycentricDynamics translated(dynamics, *source);
	const auto report = m_engine->integrateStateWithStmAndSens(translated, ng, translated.barycentric(initial, epoch),
	                                                           epoch, targets, outputs, translated.options(options));
	return translated.finish(report, initial, epoch, targets, outputs);
}
