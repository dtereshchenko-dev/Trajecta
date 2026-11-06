#pragma once
#include "IIntegrator.h"
#include "RadauCore.h"

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
