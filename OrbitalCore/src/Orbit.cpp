#include "orbital_core/Orbit.h"
#include "orbital_core/CometaryCoordinates.h"
#include <string>
#include "orbital_core/TimeScale.h"
#include <iostream>
#include <iomanip>

namespace {
	bool populateCometaryElements(Orbit& orbit)
	{
		using namespace orbital_core;
		const auto r = orbit.getR();
		const auto v = orbit.getV() * Constants::GAUSS_K;
		cometary::Vector6 state;
		state << r.x, r.y, r.z, v.x, v.y, v.z;
		const auto elements = cometary::fromCartesian(state);
		if (!elements)
			return false;

		const auto& cot = *elements;
		const double q = cot(0);
		const double e = cot(1);
		const double inverseA = (1.0 - e) / q;
		const double a = inverseA == 0.0 ? std::numeric_limits<double>::infinity() : 1.0 / inverseA;
		const auto type = e < 1.0 ? OrbitType::Elliptical : e > 1.0 ? OrbitType::Hyperbolic : OrbitType::Parabolic;
		const double nu = type == OrbitType::Elliptical ? Helpers::wrapTwoPi(cot(5)) : Helpers::wrapPi(cot(5));
		const double undefined = std::numeric_limits<double>::quiet_NaN();
		const Vec3D eccVector(e, 0.0, 0.0);
		const double anomaly =
		    type == OrbitType::Parabolic ? undefined : Calculations::excentricAnomaly(eccVector, nu, type);
		orbit.setType(type);
		orbit.setExcentricity(e);
		orbit.setPerigee(q);
		orbit.setSemiMajorAxis(a);
		orbit.setSemiMinorAxis(type == OrbitType::Parabolic ? undefined : Calculations::semiMinorAxis(a, e));
		orbit.setInclination(cot(2));
		orbit.setLongitudeOfAscendingNode(cot(3));
		orbit.setArgumentOfPericentre(cot(4));
		orbit.setTrueAnomaly(nu);
		orbit.setExcentricAnomaly(anomaly);
		orbit.setMeanAnomaly(type == OrbitType::Parabolic ? undefined
		                                                  : Calculations::meanAnomaly(eccVector, anomaly, type));
		const double absInverseA = std::abs(inverseA);
		const double meanMotion = Constants::GAUSS_K * absInverseA * std::sqrt(absInverseA);
		orbit.setMeanMotion(meanMotion);
		orbit.setOrbitalPeriod(type == OrbitType::Elliptical ? Constants::TWOPI / meanMotion : undefined);
		orbit.setApogee(type == OrbitType::Elliptical ? a * (1.0 + e) : undefined);
		orbit.setReciprocalSemiMajorAxis(inverseA);
		orbit.setTimeOfPerihelion(orbit.getEpoch() - cometary::timeSincePerihelion(q, e, cot(5)));
		return std::isfinite(orbit.getTimeOfPerihelion());
	}
}  // namespace

std::shared_ptr<Orbit> make_orbit_from_states_vectors(StateVectors vectors, const double& epoch)
{
	assert(epoch > 0);
	const Vec3D& R = vectors.first;
	const Vec3D& V = vectors.second;
	if (!std::isfinite(R.x) || !std::isfinite(R.y) || !std::isfinite(R.z) || !std::isfinite(V.x) ||
	    !std::isfinite(V.y) || !std::isfinite(V.z) || !(R.magnitude() > 0.0))
		return nullptr;
	const Vec3D& angMom = Calculations::angMomentum(R, V);
	const Vec3D& eccVec = Calculations::eccentricityVec(R, V);

	auto orbit = std::make_shared<Orbit>();
	orbit->setEpoch(epoch);
	orbit->setR(R);
	orbit->setV(V);
	if (orbital_core::cometary::usesCometaryElements(eccVec.magnitude()))
		return populateCometaryElements(*orbit) ? orbit : nullptr;

	const double& E = Calculations::specificEnergy(R, V, Constants::MU_IN_UNIT);
	const double& a = Calculations::semiMajorAxis(Constants::MU_IN_UNIT, E);
	OrbitType type = OrbitType::Undef;
	if (a > 0 && eccVec.magnitude() > 0 && eccVec.magnitude() < 1)
		type = OrbitType::Elliptical;
	else if (a < 0 && eccVec.magnitude() > 1)
		type = OrbitType::Hyperbolic;
	if (type == OrbitType::Undef)
		return nullptr;

	orbit->setType(type);
	orbit->setExcentricity(eccVec.magnitude());
	orbit->setSemiMajorAxis(a);
	orbit->setSemiMinorAxis(Calculations::semiMinorAxis(a, eccVec.magnitude()));
	orbit->setPerigee(Calculations::perigee(eccVec, orbit->getSemiMajorAxis()));
	orbit->setInclination(Calculations::inclination(angMom));
	orbit->setLongitudeOfAscendingNode(Calculations::longitudeOfAscendingNode(angMom));
	orbit->setArgumentOfPericentre(Calculations::argumentOfPericentre(eccVec, angMom));
	orbit->setTrueAnomaly(Calculations::trueAnomaly(eccVec, R, V, type));
	orbit->setExcentricAnomaly(Calculations::excentricAnomaly(eccVec, orbit->getTrueAnomaly(), type));
	orbit->setMeanAnomaly(Calculations::meanAnomaly(eccVec, orbit->getExcentricAnomaly(), type));
	orbit->setMeanMotion(Calculations::meanMotion(a, Constants::GM_SUN));
	orbit->setOrbitalPeriod(Calculations::orbitalPeriod(a, Constants::GM_SUN));
	orbit->setApogee(Calculations::apogee(eccVec, a));
	orbit->setTimeOfPerihelion(
	    Calculations::timeOfPerihelion(epoch, R, V, orbit->getMeanAnomaly(), a, Constants::GM_SUN, type));
	orbit->setReciprocalSemiMajorAxis(Calculations::reciprocalSemiMajorAxis(a));
	return orbit;
}

OrbitType Orbit::getType() const noexcept { return m_type; }

double Orbit::getEpoch() const noexcept { return m_epochMjd == 0.0 ? 0.0 : orbital_core::time::jdFromMjd(m_epochMjd); }

double Orbit::epochMjd() const noexcept { return m_epochMjd; }

double Orbit::getExcentricity() const noexcept { return m_excentricity; }

double Orbit::getSemiMajorAxis() const noexcept { return m_semiMajorAxis; }

double Orbit::getSemiMinorAxis() const noexcept { return m_semiMinorAxis; }

double Orbit::getPerigee() const noexcept { return m_perigee; }

double Orbit::getInclination() const noexcept { return m_inclination; }

double Orbit::getLongitudeOfAscendingNode() const noexcept { return m_longitudeOfAscendingNode; }

double Orbit::getArgumentOfPericentre() const noexcept { return m_argumentOfPericentre; }

double Orbit::getTrueAnomaly() const noexcept { return m_trueAnomaly; }

double Orbit::getExcentricAnomaly() const noexcept { return m_excentricAnomaly; }

double Orbit::getMeanAnomaly() const noexcept { return m_meanAnomaly; }

double Orbit::getMeanMotion() const noexcept { return m_meanMotion; }

double Orbit::getOrbitalPeriod() const noexcept { return m_orbitalPeriod; }

double Orbit::getApogee() const noexcept { return m_apogee; }

double Orbit::getTimeOfPerihelion() const noexcept { return m_timeOfPerihelion; }

double Orbit::getReciprocalSemiMajorAxis() const noexcept { return m_reciprocalSemiMajorAxis; }

Vec3D Orbit::getR() const noexcept { return m_R; }

Vec3D Orbit::getV() const noexcept { return m_V; }

std::string Orbit::getTypeString() const noexcept
{
	switch (m_type) {
		case OrbitType::Elliptical:
			return "Elliptical";
		case OrbitType::Hyperbolic:
			return "Hyperbolic";
		case OrbitType::Parabolic:
			return "Parabolic";
		default:
			return "Undefined";
	}
}

orbital_core::OrbitalClass orbital_core::orbitalClassification(const Orbit& orbit,
                                                               std::string_view designation) noexcept
{
	using Class = OrbitalClass;
	char id[32]{};
	std::size_t size = 0;
	for (const char c : designation) {
		if (c == ' ' || c == '\t')
			continue;
		if (size == sizeof(id))
			return Class::Unknown;
		id[size++] = c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : c;
	}
	std::size_t digits = 0;
	while (digits < size && id[digits] >= '0' && id[digits] <= '9')
		++digits;
	if (digits > 0 && digits < size && id[digits] == 'I' &&
	    (digits + 1 == size || id[digits + 1] == '/' || id[digits + 1] == '-'))
		return Class::Interstellar;
	const auto cometLetter = [](char c) { return c == 'P' || c == 'C' || c == 'D' || c == 'X'; };
	const bool comet = (digits > 0 && digits < size && (id[digits] == 'P' || id[digits] == 'D') &&
	                    (digits + 1 == size || id[digits + 1] == '/' || id[digits + 1] == '-' ||
	                     (digits == 4 && digits + 1 < size))) ||
	                   (size >= 2 && cometLetter(id[0]) && id[1] == '/') ||
	                   (size == 8 && cometLetter(id[0]) && (id[1] == 'I' || id[1] == 'J' || id[1] == 'K') &&
	                    id[2] >= '0' && id[2] <= '9' && id[3] >= '0' && id[3] <= '9');
	if (orbit.getType() == OrbitType::Hyperbolic && orbit.getExcentricity() > 1.0)
		return comet ? Class::HyperbolicComet : Class::Hyperbolic;
	if (orbit.getType() == OrbitType::Parabolic && orbit.getExcentricity() == 1.0)
		return comet ? Class::ParabolicComet : Class::Parabolic;
	if (comet)
		return Class::Comet;

	const double a = orbit.getSemiMajorAxis();
	const double e = orbit.getExcentricity();
	if (orbit.getType() != OrbitType::Elliptical || !std::isfinite(a) || !std::isfinite(e) || !(a > 0.0) || e < 0.0 ||
	    e >= 1.0)
		return Class::Unknown;
	const double q = a * (1.0 - e);
	const double Q = a * (1.0 + e);
	if (!std::isfinite(q) || !std::isfinite(Q))
		return Class::Unknown;
	if (a < 1.0 && Q < 0.983)
		return Class::Atira;
	if (a < 1.0 && Q >= 0.983)
		return Class::Aten;
	if (a >= 1.0 && q < 1.017)
		return Class::Apollo;
	if (a >= 1.0 && q >= 1.017 && q < 1.3)
		return Class::Amor;
	if (a < 3.2 && q >= 1.3 && q < 1.666)
		return Class::MarsCrossing;
	if (a < 2.0 && q >= 1.666)
		return Class::InnerMainBelt;
	if (a >= 2.0 && a < 3.2 && q >= 1.666)
		return Class::MainBelt;
	if (a >= 3.2 && a < 4.6)
		return Class::OuterMainBelt;
	if (a >= 4.6 && a < 5.5 && e < 0.3)
		return Class::JupiterTrojanCandidate;
	if (a >= 5.5 && a < 30.1)
		return Class::Centaur;
	if (a >= 30.1)
		return Class::Tno;
	return Class::Unknown;
}

orbital_core::CometOrbitClasses orbital_core::cometOrbitClasses(const Orbit& orbit,
                                                                std::string_view designation) noexcept
{
	CometOrbitClasses result;
	if (orbitalClassification(orbit, designation) != OrbitalClass::Comet || orbit.getType() != OrbitType::Elliptical)
		return result;
	const double a = orbit.getSemiMajorAxis();
	const double e = orbit.getExcentricity();
	if (!std::isfinite(a) || !(a > 0.0) || !std::isfinite(e) || e < 0.0 || e >= 1.0)
		return result;
	const double period = orbit.getOrbitalPeriod();
	if (std::isfinite(period) && period > 0.0) {
		if (period < 20.0 * 365.25)
			result.period = CometPeriodClass::JupiterFamily;
		else if (period <= 200.0 * 365.25)
			result.period = CometPeriodClass::HalleyType;
	}
	const double inclination = orbit.getInclination();
	if (!std::isfinite(inclination) || inclination < 0.0 || inclination > Constants::PI)
		return result;
	constexpr double jupiterAxisAu = 5.2044;
	const double tisserand =
	    jupiterAxisAu / a + 2.0 * std::cos(inclination) * std::sqrt((a / jupiterAxisAu) * (1.0 - e * e));
	if (!std::isfinite(tisserand))
		return result;
	if (tisserand > 3.0 && a < jupiterAxisAu)
		result.tisserand = CometTisserandClass::EnckeType;
	else if (tisserand > 3.0 && a > jupiterAxisAu)
		result.tisserand = CometTisserandClass::ChironType;
	else if (tisserand > 2.0 && tisserand < 3.0)
		result.tisserand = CometTisserandClass::JupiterFamily;
	return result;
}

void Orbit::setType(const OrbitType& type) noexcept { m_type = type; }

void Orbit::setEpoch(const double& epoch) noexcept
{
	m_epochMjd = epoch == 0.0 ? 0.0 : orbital_core::time::mjdFromJd(epoch);
}

void Orbit::setEpochMjd(double epochMjd) noexcept { m_epochMjd = epochMjd; }

void Orbit::setExcentricity(const double& excentricity) noexcept { m_excentricity = excentricity; }

void Orbit::setSemiMajorAxis(const double& semiMajorAxis) noexcept { m_semiMajorAxis = semiMajorAxis; }

void Orbit::setSemiMinorAxis(const double& semiMinorAxis) noexcept { m_semiMinorAxis = semiMinorAxis; }

void Orbit::setPerigee(const double& perigee) noexcept { m_perigee = perigee; }

void Orbit::setInclination(const double& inclination) noexcept { m_inclination = inclination; }

void Orbit::setLongitudeOfAscendingNode(const double& longitudeOfAscendingNode) noexcept
{
	m_longitudeOfAscendingNode = longitudeOfAscendingNode;
}

void Orbit::setArgumentOfPericentre(const double& argumentOfPericentre) noexcept
{
	m_argumentOfPericentre = argumentOfPericentre;
}

void Orbit::setTrueAnomaly(const double& trueAnomaly) noexcept { m_trueAnomaly = trueAnomaly; }

void Orbit::setExcentricAnomaly(const double& excentricAnomaly) noexcept { m_excentricAnomaly = excentricAnomaly; }

void Orbit::setMeanAnomaly(const double& meanAnomaly) noexcept { m_meanAnomaly = meanAnomaly; }

void Orbit::setMeanMotion(const double& n) noexcept { m_meanMotion = n; }

void Orbit::setOrbitalPeriod(const double& period) noexcept { m_orbitalPeriod = period; }

void Orbit::setApogee(const double& Q) noexcept { m_apogee = Q; }

void Orbit::setTimeOfPerihelion(const double& tp) noexcept { m_timeOfPerihelion = tp; }

void Orbit::setReciprocalSemiMajorAxis(const double& invA) noexcept { m_reciprocalSemiMajorAxis = invA; }

void Orbit::setR(const Vec3D& R) noexcept { m_R = R; }

void Orbit::setV(const Vec3D& V) noexcept { m_V = V; }

void Orbit::print() const noexcept { print(std::cout); }

void Orbit::print(std::ostream& os) const noexcept
{
	os << std::fixed << std::setprecision(8);
	os << "Orbit type: " << getTypeString() << "\n";
	os << "Excentricity: (e) " << m_excentricity << " \u00b1 " << getSigmaExcentricity() << "\n";
	os << "Semi-major axis: (au) " << m_semiMajorAxis << " \u00b1 " << getSigmaSemiMajorAxis() << "\n";
	os << "Perigee: (au) " << m_perigee << "\n";
	os << "Inclination: (deg) " << Convertors::radiansToDegrees(m_inclination) << " \u00b1 "
	   << Convertors::radiansToDegrees(getSigmaInclination()) << "\n";
	os << "Longitude of ascending node: (deg) " << Convertors::radiansToDegrees(m_longitudeOfAscendingNode)
	   << " \u00b1 " << Convertors::radiansToDegrees(getSigmaLongitudeOfAscendingNode()) << "\n";
	os << "Argument of pericentre: (deg) " << Convertors::radiansToDegrees(m_argumentOfPericentre) << " \u00b1 "
	   << Convertors::radiansToDegrees(getSigmaArgumentOfPericentre()) << "\n";
	os << "True anomaly: (deg) " << Convertors::radiansToDegrees(m_trueAnomaly) << " \u00b1 "
	   << Convertors::radiansToDegrees(getSigmaTrueAnomaly()) << "\n";
	os << "Excentric anomaly: (deg) " << Convertors::radiansToDegrees(m_excentricAnomaly) << "\n";
	os << "Mean anomaly: (deg) " << Convertors::radiansToDegrees(m_meanAnomaly) << "\n";
	os << "R: " << m_R << "\n";
	os << "V: " << m_V << "\n";
	os.unsetf(std::ios::fixed);
}

void Orbit::setSigmaExcentricity(const double& sigma) noexcept { m_sigmaExcentricity = sigma; }

void Orbit::setSigmaSemiMajorAxis(const double& sigma) noexcept { m_sigmaSemiMajorAxis = sigma; }

void Orbit::setSigmaInclination(const double& sigma) noexcept { m_sigmaInclination = sigma; }

void Orbit::setSigmaLongitudeOfAscendingNode(const double& sigma) noexcept { m_sigmaLongitudeOfAscendingNode = sigma; }

void Orbit::setSigmaArgumentOfPericentre(const double& sigma) noexcept { m_sigmaArgumentOfPericentre = sigma; }

void Orbit::setSigmaTrueAnomaly(const double& sigma) noexcept { m_sigmaTrueAnomaly = sigma; }

void Orbit::setSigmaPerigee(const double& sigma) noexcept { m_sigmaPerigee = sigma; }

double Orbit::getSigmaExcentricity() const noexcept { return m_sigmaExcentricity; }

double Orbit::getSigmaSemiMajorAxis() const noexcept { return m_sigmaSemiMajorAxis; }

double Orbit::getSigmaInclination() const noexcept { return m_sigmaInclination; }

double Orbit::getSigmaLongitudeOfAscendingNode() const noexcept { return m_sigmaLongitudeOfAscendingNode; }

double Orbit::getSigmaArgumentOfPericentre() const noexcept { return m_sigmaArgumentOfPericentre; }

double Orbit::getSigmaTrueAnomaly() const noexcept { return m_sigmaTrueAnomaly; }

double Orbit::getSigmaPerigee() const noexcept { return m_sigmaPerigee; }
