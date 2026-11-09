#pragma once
#include "IIntegrator.h"
#include "RadauCore.h"

namespace orbital_core::integrators::radau {
	template <int N>
	void predictRa15State(const RadauVector<N>& x, const RadauVector<N>& v, const RadauVector<N>& positionCompensation,
	                      const RadauVector<N>& velocityCompensation, const std::array<RadauVector<N>, 8>& b,
	                      const RadauVector<N>& acceleration, double fraction, double step, RadauVector<N>& outputX,
	                      RadauVector<N>& outputV)
	{
		const auto& rc = radauConstants();
		const double step2 = step * step;
		for (int k = 0; k < N; ++k) {
			double temp =
			    rc.w[3] * b[3](k) +
			    fraction *
			        (rc.w[4] * b[4](k) +
			         fraction * (rc.w[5] * b[5](k) + fraction * (rc.w[6] * b[6](k) + fraction * rc.w[7] * b[7](k))));
			outputX(k) =
			    x(k) + (fraction * (step * v(k) + step2 * fraction *
			                                          (acceleration(k) * rc.w1 +
			                                           fraction * (rc.w[1] * b[1](k) +
			                                                       fraction * (rc.w[2] * b[2](k) + fraction * temp)))) -
			            positionCompensation(k) - fraction * step * velocityCompensation(k));
			temp = rc.u[3] * b[3](k) +
			       fraction *
			           (rc.u[4] * b[4](k) +
			            fraction * (rc.u[5] * b[5](k) + fraction * (rc.u[6] * b[6](k) + fraction * rc.u[7] * b[7](k))));
			outputV(k) = v(k) + (fraction * step *
			                         (acceleration(k) + fraction * (rc.u[1] * b[1](k) +
			                                                        fraction * (rc.u[2] * b[2](k) + fraction * temp))) -
			                     velocityCompensation(k));
		}
	}

}  // namespace orbital_core::integrators::radau

// A synchronous view of a physical RA15 polynomial. Capture its low parts before
// the accepted update, then borrow coefficients only until the observer returns.
// Keeping this adapter outside the integration kernel leaves state advancement,
// requested outputs and variational propagation independent of dense consumers.
class RA15DenseOutput {
	using Vector = orbital_core::integrators::radau::RadauVector<3>;
	const AcceptedStateStepObserver& m_observer;
	const Vector* m_x = nullptr;
	const Vector* m_v = nullptr;
	const Vector* m_acceleration = nullptr;
	const std::array<Vector, 8>* m_coefficients = nullptr;
	Vector m_positionCompensation;
	Vector m_velocityCompensation;
	double m_epoch = 0, m_offset = 0, m_timeCompensation = 0, m_step = 0;
	IntegratedState interpolate(double mjd) const;

   public:
	explicit RA15DenseOutput(const AcceptedStateStepObserver& observer);
	void prepare(double epoch, double offset, double timeCompensation, double step, const Vector& x, const Vector& v,
	             const Vector& positionCompensation, const Vector& velocityCompensation, const Vector& acceleration,
	             const std::array<Vector, 8>& coefficients);
	bool operator()(double beginMjd, const Vector& beginX, const Vector& beginV, double endMjd, const Vector& endX,
	                const Vector& endV) const;
};
