#include "RA15DenseOutput.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace orbital_core::integrators::radau;

RA15DenseOutput::RA15DenseOutput(const AcceptedStateStepObserver& observer) : m_observer(observer) {}

void RA15DenseOutput::prepare(double epoch, double offset, double timeCompensation, double step, const Vector& x,
                              const Vector& v, const Vector& positionCompensation, const Vector& velocityCompensation,
                              const Vector& acceleration, const std::array<Vector, 8>& coefficients)
{
	m_epoch = epoch;
	m_offset = offset;
	m_timeCompensation = timeCompensation;
	m_step = step;
	m_x = &x;
	m_v = &v;
	m_positionCompensation = positionCompensation;
	m_velocityCompensation = velocityCompensation;
	m_acceleration = &acceleration;
	m_coefficients = &coefficients;
}

IntegratedState RA15DenseOutput::interpolate(double mjd) const
{
	const auto& rc = radauConstants();
	const auto& b = *m_coefficients;
	const auto& f1 = *m_acceleration;
	const double step = m_step, step2 = step * step;
	const double delta = ((mjd - m_epoch) - m_offset) + m_timeCompensation;
	const double s = std::clamp(delta / step, 0.0, 1.0);
	Vector x, v;
	for (int k = 0; k < 3; ++k) {
		double temp =
		    rc.w[3] * b[3](k) +
		    s * (rc.w[4] * b[4](k) + s * (rc.w[5] * b[5](k) + s * (rc.w[6] * b[6](k) + s * rc.w[7] * b[7](k))));
		x(k) = (*m_x)(k) +
		       (s * (step * (*m_v)(k) +
		             step2 * s * (f1(k) * rc.w1 + s * (rc.w[1] * b[1](k) + s * (rc.w[2] * b[2](k) + s * temp)))) -
		        m_positionCompensation(k) - s * step * m_velocityCompensation(k));
		temp = rc.u[3] * b[3](k) +
		       s * (rc.u[4] * b[4](k) + s * (rc.u[5] * b[5](k) + s * (rc.u[6] * b[6](k) + s * rc.u[7] * b[7](k))));
		v(k) = (*m_v)(k) + (s * step * (f1(k) + s * (rc.u[1] * b[1](k) + s * (rc.u[2] * b[2](k) + s * temp))) -
		                    m_velocityCompensation(k));
	}
	return {positionFromRadau(x), velocityFromRadau(v)};
}

bool RA15DenseOutput::operator()(double beginMjd, const Vector& beginX, const Vector& beginV, double endMjd,
                                 const Vector& endX, const Vector& endV) const
{
	AcceptedStateStep step{beginMjd,
	                       endMjd,
	                       {positionFromRadau(beginX), velocityFromRadau(beginV)},
	                       {positionFromRadau(endX), velocityFromRadau(endV)}};
	step.interpolateMjd = [this, &step](double mjd) {
		if (!std::isfinite(mjd) || mjd < std::min(step.beginMjd, step.endMjd) ||
		    mjd > std::max(step.beginMjd, step.endMjd))
			throw std::out_of_range("RA15 dense output is outside the accepted step");
		if (mjd == step.beginMjd)
			return step.begin;
		if (mjd == step.endMjd)
			return step.end;
		return interpolate(mjd);
	};
	return m_observer(step);
}
