#ifndef ORBIT_H
#define ORBIT_H
#include "orbital_core/Vec3D.h"
#include <string>
#include <string_view>
#include "orbital_core/Types.h"
#include "orbital_core/StateVectors.h"
#include <cassert>
#include <memory>

class ORBITAL_CORE_API Orbit {
   public:
	Orbit() = default;
	Orbit(const Orbit& other) = default;
	~Orbit() = default;

	[[nodiscard]] OrbitType getType() const noexcept;
	[[nodiscard]] double getEpoch() const noexcept;
	[[nodiscard]] double getExcentricity() const noexcept;
	[[nodiscard]] double getSemiMajorAxis() const noexcept;
	[[nodiscard]] double getSemiMinorAxis() const noexcept;
	[[nodiscard]] double getPerigee() const noexcept;
	[[nodiscard]] double getInclination() const noexcept;
	[[nodiscard]] double getLongitudeOfAscendingNode() const noexcept;
	[[nodiscard]] double getArgumentOfPericentre() const noexcept;
	[[nodiscard]] double getTrueAnomaly() const noexcept;
	[[nodiscard]] double getExcentricAnomaly() const noexcept;
	[[nodiscard]] double getMeanAnomaly() const noexcept;
	[[nodiscard]] double getMeanMotion() const noexcept;               // n = sqrt(µ/|a|³)
	[[nodiscard]] double getOrbitalPeriod() const noexcept;            // P = 2π/n (days)
	[[nodiscard]] double getApogee() const noexcept;                   // Q = a(1+e)
	[[nodiscard]] double getTimeOfPerihelion() const noexcept;         // Tp (days, same epoch frame)
	[[nodiscard]] double getReciprocalSemiMajorAxis() const noexcept;  // 1/a
	[[nodiscard]] std::string getTypeString() const noexcept;

	[[nodiscard]] Vec3D getR() const noexcept;
	[[nodiscard]] Vec3D getV() const noexcept;

	void setType(const OrbitType& type) noexcept;
	// Epoch. Stored as MJD TDB: a double holds MJD 59169 with an ulp of
	// 0.63 us against 40.23 us at the matching JD, and only 1 of every 64
	// representable instants survives a round trip through JD. getEpoch()
	// keeps returning the JD view so existing readers are unaffected; reach
	// for epochMjd() wherever the resolution matters. See TimeScale.h.
	void setEpoch(const double& epoch) noexcept;  // takes JD
	void setEpochMjd(double epochMjd) noexcept;   // lossless
	[[nodiscard]] double epochMjd() const noexcept;
	void setExcentricity(const double& excentricity) noexcept;
	void setSemiMajorAxis(const double& semiMajorAxis) noexcept;
	void setSemiMinorAxis(const double& semiMinorAxis) noexcept;
	void setPerigee(const double& perigee) noexcept;
	void setInclination(const double& inclination) noexcept;
	void setLongitudeOfAscendingNode(const double& longitudeOfAscendingNode) noexcept;
	void setArgumentOfPericentre(const double& argumentOfPericentre) noexcept;
	void setTrueAnomaly(const double& trueAnomaly) noexcept;
	void setExcentricAnomaly(const double& excentricAnomaly) noexcept;
	void setMeanAnomaly(const double& meanAnomaly) noexcept;
	void setMeanMotion(const double& n) noexcept;
	void setOrbitalPeriod(const double& period) noexcept;
	void setApogee(const double& Q) noexcept;
	void setTimeOfPerihelion(const double& tp) noexcept;
	void setReciprocalSemiMajorAxis(const double& invA) noexcept;

	// NEW
	void setR(const Vec3D& R) noexcept;
	void setV(const Vec3D& V) noexcept;

	void print() const noexcept;
	void print(std::ostream& os) const noexcept;

	void setSigmaExcentricity(const double& sigma) noexcept;
	void setSigmaSemiMajorAxis(const double& sigma) noexcept;
	void setSigmaInclination(const double& sigma) noexcept;
	void setSigmaLongitudeOfAscendingNode(const double& sigma) noexcept;
	void setSigmaArgumentOfPericentre(const double& sigma) noexcept;
	void setSigmaTrueAnomaly(const double& sigma) noexcept;
	void setSigmaPerigee(const double& sigma) noexcept;

	[[nodiscard]] double getSigmaExcentricity() const noexcept;
	[[nodiscard]] double getSigmaSemiMajorAxis() const noexcept;
	[[nodiscard]] double getSigmaInclination() const noexcept;
	[[nodiscard]] double getSigmaLongitudeOfAscendingNode() const noexcept;
	[[nodiscard]] double getSigmaArgumentOfPericentre() const noexcept;
	[[nodiscard]] double getSigmaTrueAnomaly() const noexcept;
	[[nodiscard]] double getSigmaPerigee() const noexcept;

   private:
	OrbitType m_type = OrbitType::Undef;
	double m_epochMjd{};
	double m_excentricity{};
	double m_semiMajorAxis{};
	double m_semiMinorAxis{};
	double m_perigee{};
	double m_inclination{};
	double m_longitudeOfAscendingNode{};
	double m_argumentOfPericentre{};
	double m_trueAnomaly{};
	double m_excentricAnomaly{};
	double m_meanAnomaly{};
	double m_meanMotion{};
	double m_orbitalPeriod{};
	double m_apogee{};
	double m_timeOfPerihelion{};
	double m_reciprocalSemiMajorAxis{};

	double m_sigmaExcentricity{};
	double m_sigmaSemiMajorAxis{};
	double m_sigmaInclination{};
	double m_sigmaLongitudeOfAscendingNode{};
	double m_sigmaArgumentOfPericentre{};
	double m_sigmaTrueAnomaly{};
	double m_sigmaPerigee{};

	// NEW
	Vec3D m_R{};
	Vec3D m_V{};
};

[[nodiscard]] ORBITAL_CORE_API std::shared_ptr<Orbit> make_orbit_from_states_vectors(StateVectors vectors,
                                                                                     const double& epoch);

using Orbit_p = std::shared_ptr<Orbit>;

namespace orbital_core {
	enum class OrbitalClass {
		Unknown,
		Interstellar,
		Comet,
		HyperbolicComet,
		ParabolicComet,
		Hyperbolic,
		Parabolic,
		Atira,
		Aten,
		Apollo,
		Amor,
		MarsCrossing,
		InnerMainBelt,
		MainBelt,
		OuterMainBelt,
		JupiterTrojanCandidate,
		Centaur,
		Tno,
	};

	enum class CometPeriodClass { Unknown, JupiterFamily, HalleyType };
	enum class CometTisserandClass { Unknown, EnckeType, ChironType, JupiterFamily };

	struct CometOrbitClasses {
		CometPeriodClass period{CometPeriodClass::Unknown};
		CometTisserandClass tisserand{CometTisserandClass::Unknown};
	};

	// An orbital group is not a physical identification; explicit comet and
	// interstellar designations take precedence over the NEO orbital groups.
	[[nodiscard]] ORBITAL_CORE_API OrbitalClass orbitalClassification(const Orbit& orbit,
	                                                                  std::string_view designation) noexcept;
	[[nodiscard]] ORBITAL_CORE_API CometOrbitClasses cometOrbitClasses(const Orbit& orbit,
	                                                                   std::string_view designation) noexcept;
}  // namespace orbital_core

#endif  // ORBIT_H
