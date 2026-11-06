#include "orbital_core/nongrav/NonGravModel.h"
#include <cstdint>
#include "orbital_core/Types.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>

namespace {

	struct MarsdenCoefficients {
		double alpha;
		double r0;
		double m;
		double n;
		double k;
	};

	inline constexpr MarsdenCoefficients kMarsdenH2O{0.1112620426, 2.808, 2.15, 5.093, 4.6142};

	inline constexpr MarsdenCoefficients kMarsdenCO{0.0408373333128795, 5.0, 2.0, 3.0, 2.6};

	[[nodiscard]] double evalMarsden(double rAU, const MarsdenCoefficients& c) noexcept
	{
		const double rr = rAU / c.r0;
		return c.alpha * std::pow(rr, -c.m) * std::pow(1.0 + std::pow(rr, c.n), -c.k);
	}

}  // namespace

void NonGravModel::setModel(NgModel m) noexcept
{
	m_model = m;
	m_fitMask = defaultFitParameterMaskFor(m);
}

NgParameterMask NonGravModel::defaultFitParameterMaskFor(NgModel model) noexcept
{
	switch (model) {
		case NgModel::None:
			return NgParameterMask::None;

		case NgModel::Yarkovsky:
			return NgParameterMask::A2;

		case NgModel::InverseSquare:
			return NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3;

		case NgModel::MarsdenH2O:
		case NgModel::MarsdenCO:
			return NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 | NgParameterMask::DT;
	}

	return NgParameterMask::None;
}

NgParameterMask NonGravModel::sanitizeFitParameterMaskFor(NgModel model, NgParameterMask mask) noexcept
{
	if (mask == NgParameterMask::None)
		return NgParameterMask::None;

	constexpr NgParameterMask kA1A2 = NgParameterMask::A1 | NgParameterMask::A2;
	constexpr NgParameterMask kA1A2A3 = NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3;
	constexpr NgParameterMask kA1A2A3DT =
	    NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 | NgParameterMask::DT;

	switch (model) {
		case NgModel::None:
			return NgParameterMask::None;

		case NgModel::Yarkovsky:
			return hasNgParameter(mask, NgParameterMask::A2) ? NgParameterMask::A2 : NgParameterMask::None;
		case NgModel::InverseSquare:
		case NgModel::MarsdenH2O:
		case NgModel::MarsdenCO:
			if (hasNgParameter(mask, NgParameterMask::DT))
				return kA1A2A3DT;

			if (hasNgParameter(mask, NgParameterMask::A3))
				return kA1A2A3;

			if (hasNgParameter(mask, NgParameterMask::A1) || hasNgParameter(mask, NgParameterMask::A2))
				return kA1A2;

			return NgParameterMask::None;
	}

	return NgParameterMask::None;
}

void NonGravModel::setFitParameterMask(NgParameterMask mask) noexcept
{
	m_fitMask = sanitizeFitParameterMaskFor(m_model, mask);
}

NgParameterMask NonGravModel::fitParameterMask() const noexcept
{
	if (!enabled())
		return NgParameterMask::None;

	return sanitizeFitParameterMaskFor(m_model, m_fitMask);
}

int NonGravModel::fitParameterCount() const noexcept
{
	const std::uint8_t m = static_cast<std::uint8_t>(fitParameterMask());

	int count = 0;
	for (int bit = 0; bit < maxFitParams; ++bit)
		count += (m >> bit) & 1u;

	return count;
}

NgParameter NonGravModel::fitParameterAt(int fitIndex) const noexcept
{
	const NgParameterMask mask = fitParameterMask();
	int seen = 0;
	const NgParameter order[maxFitParams] = {NgParameter::A1, NgParameter::A2, NgParameter::A3, NgParameter::DT};

	const NgParameterMask bits[maxFitParams] = {NgParameterMask::A1, NgParameterMask::A2, NgParameterMask::A3,
	                                            NgParameterMask::DT};

	for (int i = 0; i < maxFitParams; ++i) {
		if (!hasNgParameter(mask, bits[i]))
			continue;

		if (seen == fitIndex)
			return order[i];

		++seen;
	}
	return NgParameter::A1;
}

int NonGravModel::sensitivityColumnForFitParameter(int fitIndex) const noexcept
{
	switch (fitParameterAt(fitIndex)) {
		case NgParameter::A1:
			return 0;
		case NgParameter::A2:
			return 1;
		case NgParameter::A3:
			return 2;
		case NgParameter::DT:
			return 3;
	}
	return 0;
}

double NonGravModel::fitParameterValue(int fitIndex) const noexcept
{
	switch (fitParameterAt(fitIndex)) {
		case NgParameter::A1:
			return m_A1;
		case NgParameter::A2:
			return m_A2;
		case NgParameter::A3:
			return m_A3;
		case NgParameter::DT:
			return m_DT;
	}
	return 0.0;
}

void NonGravModel::setFitParameterValue(int fitIndex, double value) noexcept
{
	switch (fitParameterAt(fitIndex)) {
		case NgParameter::A1:
			m_A1 = value;
			break;
		case NgParameter::A2:
			m_A2 = value;
			break;
		case NgParameter::A3:
			m_A3 = value;
			break;
		case NgParameter::DT:
			m_DT = value;
			break;
	}
}

double NonGravModel::gr(double rAU) const noexcept
{
	if (rAU <= 0.0)
		return 0.0;

	switch (m_model) {
		case NgModel::None:
			return 0.0;

		case NgModel::MarsdenH2O:
			return evalMarsden(rAU, kMarsdenH2O);

		case NgModel::MarsdenCO:
			return evalMarsden(rAU, kMarsdenCO);

		case NgModel::InverseSquare:
		case NgModel::Yarkovsky:
			return 1.0 / (rAU * rAU);
	}

	return 0.0;
}

double NonGravModel::rPastMagnitude(const Vec3D& r, const Vec3D& v, double dtDays) const noexcept
{
	const double rMag = r.magnitude();
	if (dtDays == 0.0)
		return rMag;

	constexpr double mu = Constants::GM_SUN;

	const double eps = Calculations::specificEnergy(r, v, mu);
	if (std::fabs(eps) < 1e-18)
		return rMag;

	const double a = Calculations::semiMajorAxis(mu, eps);
	const double h = Calculations::angMomentum(r, v).magnitude();

	const double rv = r.dot(v);
	const OrbitType type = (eps < 0.0) ? OrbitType::Elliptical : OrbitType::Hyperbolic;

	if (type == OrbitType::Elliptical) {
		const double arg = 1.0 - h * h / (mu * a);
		const double e = std::sqrt(std::max(0.0, arg));
		if (e < 1e-12)
			return rMag;

		double cosE = std::clamp((1.0 - rMag / a) / e, -1.0, 1.0);
		double E_now = std::acos(cosE);
		if (rv < 0.0)
			E_now = -E_now;

		const double M_now = E_now - e * std::sin(E_now);
		const double n_mean = std::sqrt(mu / (a * a * a));
		const double M_past = M_now - n_mean * dtDays;

		const double E_past = Calculations::kepler_solve_e(M_past, e);
		return a * (1.0 - e * std::cos(E_past));
	}

	const double absA = std::fabs(a);
	const double arg = 1.0 + h * h / (mu * absA);
	const double e = std::sqrt(std::max(1.0 + 1e-18, arg));

	double coshF = (1.0 - rMag / a) / e;
	if (coshF < 1.0)
		coshF = 1.0;

	double F_now = std::acosh(coshF);
	if (rv < 0.0)
		F_now = -F_now;

	const double M_now = e * std::sinh(F_now) - F_now;
	const double n_mean = std::sqrt(mu / (absA * absA * absA));
	const double M_past = M_now - n_mean * dtDays;

	const double F_past = Calculations::kepler_solve_h(M_past, e);
	return a * (1.0 - e * std::cosh(F_past));
}

double NonGravModel::grDerivative(double rAU) const noexcept
{
	if (!(rAU > 0.0) || m_model == NgModel::None)
		return 0.0;
	if (m_model == NgModel::MarsdenH2O || m_model == NgModel::MarsdenCO) {
		const auto& c = m_model == NgModel::MarsdenH2O ? kMarsdenH2O : kMarsdenCO;
		const double power = std::pow(rAU / c.r0, c.n);
		return -gr(rAU) * (c.m + c.k * c.n * power / (1.0 + power)) / rAU;
	}
	return -2.0 * gr(rAU) / rAU;
}

Eigen::Matrix<double, 1, 6> NonGravModel::rPastGradient(const Vec3D& r, const Vec3D& v, double dtDays) const noexcept
{
	using Row = Eigen::Matrix<double, 1, 6>;
	const Eigen::Vector3d R(r.x, r.y, r.z), V(v.x, v.y, v.z);
	const double radius = r.magnitude();
	Row dr = Row::Zero();
	dr.head<3>() = R.transpose() / radius;
	constexpr double mu = Constants::GM_SUN;
	const double energy = Calculations::specificEnergy(r, v, mu);
	if (dtDays == 0.0 || std::abs(energy) < 1e-18)
		return dr;
	const double a = Calculations::semiMajorAxis(mu, energy);
	const double h = Calculations::angMomentum(r, v).magnitude();
	const double sigma = r.dot(v);
	Row dEnergy, dSigma, dSpeed2;
	dEnergy << (mu / (radius * radius * radius) * R).transpose(), V.transpose();
	dSigma << V.transpose(), R.transpose();
	dSpeed2 << 0.0, 0.0, 0.0, (2.0 * V).transpose();
	const Row da = mu / (2.0 * energy * energy) * dEnergy;
	const Row dh2 = 2.0 * radius * v.dot(v) * dr + radius * radius * dSpeed2 - 2.0 * sigma * dSigma;
	const bool elliptic = energy < 0.0;
	const double e = std::sqrt(std::max(elliptic ? 0.0 : 1.0 + 1e-18, 1.0 - h * h / (mu * a)));
	if (elliptic && e < 1e-12)
		return dr;
	const Row de = (-dh2 / (mu * a) + h * h / (mu * a * a) * da) / (2.0 * e);
	const double absA = std::abs(a);
	const Row dAbsA = elliptic ? da : Row(-da);
	const double root = std::sqrt(mu * absA);
	const double s = sigma / root;
	const double c = 1.0 - radius / a;
	const Row ds = dSigma / root - s / (2.0 * absA) * dAbsA;
	const Row dc = -dr / a + radius / (a * a) * da;
	const Row dAnomaly = (c * ds - s * dc) / (e * e);
	const double n = std::sqrt(mu / (absA * absA * absA));
	const Row dn = -1.5 * n / absA * dAbsA;
	if (elliptic) {
		double anomaly = std::acos(std::clamp(c / e, -1.0, 1.0));
		if (sigma < 0.0)
			anomaly = -anomaly;
		const double mean = anomaly - e * std::sin(anomaly) - n * dtDays;
		const Row dMean = (1.0 - c) * dAnomaly - std::sin(anomaly) * de - dtDays * dn;
		const double past = Calculations::kepler_solve_e(mean, e);
		const double cosPast = std::cos(past), sinPast = std::sin(past);
		const double factor = 1.0 - e * cosPast;
		const Row dPast = (dMean + sinPast * de) / factor;
		return factor * da + a * (-cosPast * de + e * sinPast * dPast);
	}
	double anomaly = std::acosh(std::max(1.0, c / e));
	if (sigma < 0.0)
		anomaly = -anomaly;
	const double mean = e * std::sinh(anomaly) - anomaly - n * dtDays;
	const Row dMean = (c - 1.0) * dAnomaly + std::sinh(anomaly) * de - dtDays * dn;
	const double past = Calculations::kepler_solve_h(mean, e);
	const double coshPast = std::cosh(past), sinhPast = std::sinh(past);
	const Row dPast = (dMean - sinhPast * de) / (e * coshPast - 1.0);
	return (1.0 - e * coshPast) * da + a * (-coshPast * de - e * sinhPast * dPast);
}

void NonGravModel::accelerationPartials(const Vec3D& r, const Vec3D& v, Eigen::Matrix3d& dadr,
                                        Eigen::Matrix3d& dadv) const noexcept
{
	dadr.setZero();
	dadv.setZero();
	const double radius = r.magnitude();
	if (!enabled() || radius < 1e-12 || (m_A1 == 0.0 && m_A2 == 0.0 && m_A3 == 0.0))
		return;
	const Eigen::Vector3d R(r.x, r.y, r.z), V(v.x, v.y, v.z);
	const Eigen::Vector3d radial = R / radius;
	const Eigen::Matrix3d dradial = (Eigen::Matrix3d::Identity() - radial * radial.transpose()) / radius;
	Eigen::Vector3d direction = m_A1 * radial;
	Eigen::Matrix3d directionR = m_A1 * dradial;
	Eigen::Matrix3d directionV = Eigen::Matrix3d::Zero();
	const auto crossMatrix = [](const Eigen::Vector3d& x) {
		Eigen::Matrix3d out;
		out << 0.0, -x.z(), x.y(), x.z(), 0.0, -x.x(), -x.y(), x.x(), 0.0;
		return out;
	};
	const auto h = r.cross(v);
	const double hMag = h.magnitude();
	if (hMag >= 1e-15) {
		const Eigen::Vector3d normal(h.x / hMag, h.y / hMag, h.z / hMag);
		const Eigen::Vector3d transverse = normal.cross(radial);
		const Eigen::Matrix3d projection = (Eigen::Matrix3d::Identity() - normal * normal.transpose()) / hMag;
		const Eigen::Matrix3d normalR = -projection * crossMatrix(V);
		const Eigen::Matrix3d normalV = projection * crossMatrix(R);
		const Eigen::Matrix3d transverseR = -crossMatrix(radial) * normalR + crossMatrix(normal) * dradial;
		const Eigen::Matrix3d transverseV = -crossMatrix(radial) * normalV;
		direction += m_A2 * transverse + m_A3 * normal;
		directionR += m_A2 * transverseR + m_A3 * normalR;
		directionV = m_A2 * transverseV + m_A3 * normalV;
	}
	const double dt = usesDT() ? m_DT : 0.0;
	const double pastRadius = dt == 0.0 ? radius : rPastMagnitude(r, v, dt);
	const auto gradient = rPastGradient(r, v, dt);
	const double g = gr(pastRadius), dg = grDerivative(pastRadius);
	dadr = g * directionR + dg * direction * gradient.head<3>();
	dadv = g * directionV + dg * direction * gradient.tail<3>();
}

double NonGravModel::rPastDerivDT(const Vec3D& r, const Vec3D& v, double dtDays) const noexcept
{
	const double rMag = r.magnitude();
	if (dtDays == 0.0) {
		return (rMag > 1e-12) ? -r.dot(v) / rMag : 0.0;
	}

	constexpr double mu = Constants::GM_SUN;
	const double eps = Calculations::specificEnergy(r, v, mu);
	if (std::fabs(eps) < 1e-18)
		return 0.0;

	const double a = Calculations::semiMajorAxis(mu, eps);
	const double hMom = Calculations::angMomentum(r, v).magnitude();
	const double rv = r.dot(v);
	const OrbitType type = (eps < 0.0) ? OrbitType::Elliptical : OrbitType::Hyperbolic;

	if (type == OrbitType::Elliptical) {
		const double arg = 1.0 - hMom * hMom / (mu * a);
		const double e = std::sqrt(std::max(0.0, arg));
		if (e < 1e-12)
			return 0.0;

		double cosE_now = std::clamp((1.0 - rMag / a) / e, -1.0, 1.0);
		double E_now = std::acos(cosE_now);
		if (rv < 0.0)
			E_now = -E_now;

		const double M_now = E_now - e * std::sin(E_now);
		const double n_mean = std::sqrt(mu / (a * a * a));
		const double M_past = M_now - n_mean * dtDays;
		const double E_past = Calculations::kepler_solve_e(M_past, e);

		const double denom = 1.0 - e * std::cos(E_past);
		if (std::fabs(denom) < 1e-15)
			return 0.0;

		return -a * e * n_mean * std::sin(E_past) / denom;
	}

	const double absA = std::fabs(a);
	const double arg = 1.0 + hMom * hMom / (mu * absA);
	const double e = std::sqrt(std::max(1.0 + 1e-18, arg));

	double coshF_now = (1.0 - rMag / a) / e;
	if (coshF_now < 1.0)
		coshF_now = 1.0;
	double F_now = std::acosh(coshF_now);
	if (rv < 0.0)
		F_now = -F_now;

	const double M_now = e * std::sinh(F_now) - F_now;
	const double n_mean = std::sqrt(mu / (absA * absA * absA));
	const double M_past = M_now - n_mean * dtDays;
	const double F_past = Calculations::kepler_solve_h(M_past, e);

	const double denom = e * std::cosh(F_past) - 1.0;
	if (std::fabs(denom) < 1e-15)
		return 0.0;

	return a * e * n_mean * std::sinh(F_past) / denom;
}

Vec3D NonGravModel::sensitivityPartialDT(const Vec3D& r, const Vec3D& v, double) const noexcept
{
	if (!enabled() || !usesDT())
		return Vec3D();

	const double rMag = r.magnitude();
	if (rMag < 1e-12)
		return Vec3D();

	const double rForG = (m_DT != 0.0) ? rPastMagnitude(r, v, m_DT) : rMag;
	const double dRForG_dDT = rPastDerivDT(r, v, m_DT);

	const double h_g = std::max(std::abs(rForG) * 1e-6, 1e-10);
	const double g_prime = (gr(rForG + h_g) - gr(rForG - h_g)) / (2.0 * h_g);

	const Vec3D rHat = r * (1.0 / rMag);
	const Vec3D hvec = r.cross(v);
	const double hMag = hvec.magnitude();
	Vec3D direction;
	if (hMag < 1e-15) {
		direction = rHat * m_A1;
	} else {
		const Vec3D nHat = hvec * (1.0 / hMag);
		const Vec3D tHat = nHat.cross(rHat);
		direction = rHat * m_A1 + tHat * m_A2 + nHat * m_A3;
	}

	return direction * (g_prime * dRForG_dDT);
}

void NonGravModel::sensitivityPartialsA(const Vec3D& r, const Vec3D& v, Vec3D& dA1, Vec3D& dA2,
                                        Vec3D& dA3) const noexcept
{
	dA1 = dA2 = dA3 = Vec3D();
	if (!enabled())
		return;

	const double rMag = r.magnitude();
	if (rMag < 1e-12)
		return;

	double rForG = rMag;
	if (usesDT() && m_DT != 0.0)
		rForG = rPastMagnitude(r, v, m_DT);
	const double g = gr(rForG);

	const Vec3D rHat = r * (1.0 / rMag);

	const Vec3D h = r.cross(v);
	const double hMag = h.magnitude();
	if (hMag < 1e-15) {
		dA1 = rHat * g;
		return;
	}

	const Vec3D nHat = h * (1.0 / hMag);
	const Vec3D tHat = nHat.cross(rHat);
	dA1 = rHat * g;
	dA2 = tHat * g;
	dA3 = nHat * g;
}

Vec3D NonGravModel::acceleration(const Vec3D& r, const Vec3D& v, double) const noexcept
{
	if (!enabled())
		return Vec3D();
	if (m_A1 == 0.0 && m_A2 == 0.0 && m_A3 == 0.0)
		return Vec3D();

	const double rMag = r.magnitude();
	if (rMag < 1e-12)
		return Vec3D();

	double rForG = rMag;
	if (usesDT() && m_DT != 0.0)
		rForG = rPastMagnitude(r, v, m_DT);

	const double g = gr(rForG);
	const Vec3D rHat = r * (1.0 / rMag);

	const Vec3D h = r.cross(v);
	const double hMag = h.magnitude();
	if (hMag < 1e-15) {
		return rHat * (m_A1 * g);
	}

	const Vec3D nHat = h * (1.0 / hMag);
	const Vec3D tHat = nHat.cross(rHat);
	return rHat * (m_A1 * g) + tHat * (m_A2 * g) + nHat * (m_A3 * g);
}

int NonGravModel::numFitParams() const noexcept { return fitParameterCount(); }
