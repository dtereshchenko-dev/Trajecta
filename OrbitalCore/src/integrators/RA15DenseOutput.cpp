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
	const double delta = ((mjd - m_epoch) - m_offset) + m_timeCompensation;
	const double s = std::clamp(delta / m_step, 0.0, 1.0);
	Vector x, v;
	predictRa15State(*m_x, *m_v, m_positionCompensation, m_velocityCompensation, *m_coefficients, *m_acceleration,
	                 s, m_step, x, v);
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
