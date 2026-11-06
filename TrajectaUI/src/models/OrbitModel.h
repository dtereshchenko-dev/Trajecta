#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QString>
#include <QStringList>

#include "orbital_core/Orbit.h"
#include "orbital_core/uncertainty/DerivedSigmas.h"

// Cached orbit values used by the workspace. Angles and their sigmas are in degrees.
class OrbitModel : public QObject {
	Q_OBJECT
	QML_ELEMENT

	Q_PROPERTY(QString orbitType READ orbitType NOTIFY changed)
	Q_PROPERTY(QString orbitalClass READ orbitalClass NOTIFY changed)
	Q_PROPERTY(QStringList orbitalClasses READ orbitalClasses NOTIFY changed)
	Q_PROPERTY(double eccentricity READ eccentricity NOTIFY changed)
	Q_PROPERTY(double semiMajorAxis READ semiMajorAxis NOTIFY changed)
	Q_PROPERTY(double perigee READ perigee NOTIFY changed)
	Q_PROPERTY(double apogee READ apogee NOTIFY changed)
	Q_PROPERTY(double inclination READ inclination NOTIFY changed)
	Q_PROPERTY(double ascendingNode READ ascendingNode NOTIFY changed)
	Q_PROPERTY(double argPericentre READ argPericentre NOTIFY changed)
	Q_PROPERTY(double trueAnomaly READ trueAnomaly NOTIFY changed)
	Q_PROPERTY(double meanAnomaly READ meanAnomaly NOTIFY changed)
	Q_PROPERTY(double meanMotion READ meanMotion NOTIFY changed)
	Q_PROPERTY(double orbitalPeriod READ orbitalPeriod NOTIFY changed)
	Q_PROPERTY(double reciprocalSemiMajorAxis READ reciprocalSemiMajorAxis NOTIFY changed)
	Q_PROPERTY(double epoch READ epoch NOTIFY changed)

	Q_PROPERTY(double sigmaE READ sigmaE NOTIFY changed)
	Q_PROPERTY(double sigmaA READ sigmaA NOTIFY changed)
	Q_PROPERTY(double sigmaI READ sigmaI NOTIFY changed)
	Q_PROPERTY(double sigmaNode READ sigmaNode NOTIFY changed)
	Q_PROPERTY(double sigmaW READ sigmaW NOTIFY changed)
	Q_PROPERTY(double sigmaNu READ sigmaNu NOTIFY changed)
	Q_PROPERTY(double sigmaM READ sigmaM NOTIFY changed)
	Q_PROPERTY(double sigmaQ READ sigmaQ NOTIFY changed)
	Q_PROPERTY(double sigmaPeriod READ sigmaPeriod NOTIFY changed)

   public:
	explicit OrbitModel(QObject* parent = nullptr);
	~OrbitModel() override;

	// Single mutator.  `derivedSigmas` may be empty (default-constructed
	// DerivedSigmas{}); pass `hasDerivedSigmas = false` to suppress the
	// sigma readouts on fits without computed covariance.
	void setFromFit(const Orbit_p& orbit, const orbital_core::uncertainty::DerivedSigmas& derivedSigmas,
	                bool hasDerivedSigmas,
	                orbital_core::OrbitalClass orbitalClass = orbital_core::OrbitalClass::Unknown,
	                orbital_core::CometOrbitClasses cometClasses = {});

	void clear();

	[[nodiscard]] QString orbitType() const noexcept { return m_orbitType; }

	[[nodiscard]] QString orbitalClass() const noexcept { return m_orbitalClass; }

	[[nodiscard]] QStringList orbitalClasses() const { return m_orbitalClasses; }

	[[nodiscard]] double eccentricity() const noexcept { return m_e; }

	[[nodiscard]] double semiMajorAxis() const noexcept { return m_a; }

	[[nodiscard]] double perigee() const noexcept { return m_q; }

	[[nodiscard]] double apogee() const noexcept { return m_Q; }

	[[nodiscard]] double inclination() const noexcept { return m_i; }

	[[nodiscard]] double ascendingNode() const noexcept { return m_node; }

	[[nodiscard]] double argPericentre() const noexcept { return m_w; }

	[[nodiscard]] double trueAnomaly() const noexcept { return m_nu; }

	[[nodiscard]] double meanAnomaly() const noexcept { return m_M; }

	[[nodiscard]] double meanMotion() const noexcept { return m_n; }

	[[nodiscard]] double orbitalPeriod() const noexcept { return m_period; }

	[[nodiscard]] double reciprocalSemiMajorAxis() const noexcept { return m_invA; }

	[[nodiscard]] double epoch() const noexcept { return m_epoch; }

	[[nodiscard]] double sigmaE() const noexcept { return m_sigmaE; }

	[[nodiscard]] double sigmaA() const noexcept { return m_sigmaA; }

	[[nodiscard]] double sigmaI() const noexcept { return m_sigmaI; }

	[[nodiscard]] double sigmaNode() const noexcept { return m_sigmaNode; }

	[[nodiscard]] double sigmaW() const noexcept { return m_sigmaW; }

	[[nodiscard]] double sigmaNu() const noexcept { return m_sigmaNu; }

	[[nodiscard]] double sigmaM() const noexcept { return m_sigmaM; }

	[[nodiscard]] double sigmaQ() const noexcept { return m_sigmaQ; }

	[[nodiscard]] double sigmaPeriod() const noexcept { return m_sigmaPeriod; }

   signals:
	void changed();

   private:
	QString m_orbitType;
	QString m_orbitalClass;
	QStringList m_orbitalClasses;
	double m_e{}, m_a{}, m_q{}, m_Q{};
	double m_i{}, m_node{}, m_w{}, m_nu{}, m_M{};
	double m_n{}, m_period{}, m_invA{}, m_epoch{};

	double m_sigmaE{}, m_sigmaA{}, m_sigmaI{}, m_sigmaNode{};
	double m_sigmaW{}, m_sigmaNu{}, m_sigmaM{}, m_sigmaQ{};
	double m_sigmaPeriod{};
};
