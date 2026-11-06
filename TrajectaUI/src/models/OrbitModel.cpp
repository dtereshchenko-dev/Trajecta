#include "OrbitModel.h"

#include "orbital_core/Types.h"

namespace {
	QString classLabel(orbital_core::OrbitalClass value)
	{
		using C = orbital_core::OrbitalClass;
		switch (value) {
			case C::Unknown:
				return {};
			case C::Interstellar:
				return OrbitModel::tr("Interstellar");
			case C::Comet:
				return OrbitModel::tr("Comet");
			case C::HyperbolicComet:
				return OrbitModel::tr("Hyperbolic comet");
			case C::ParabolicComet:
				return OrbitModel::tr("Parabolic comet");
			case C::Hyperbolic:
				return OrbitModel::tr("Hyperbolic");
			case C::Parabolic:
				return OrbitModel::tr("Parabolic");
			case C::Atira:
				return OrbitModel::tr("Atira");
			case C::Aten:
				return OrbitModel::tr("Aten");
			case C::Apollo:
				return OrbitModel::tr("Apollo");
			case C::Amor:
				return OrbitModel::tr("Amor");
			case C::MarsCrossing:
				return OrbitModel::tr("Mars-crossing");
			case C::InnerMainBelt:
				return OrbitModel::tr("Inner main belt");
			case C::MainBelt:
				return OrbitModel::tr("Main belt");
			case C::OuterMainBelt:
				return OrbitModel::tr("Outer main belt");
			case C::JupiterTrojanCandidate:
				return OrbitModel::tr("Jupiter Trojan candidate");
			case C::Centaur:
				return OrbitModel::tr("Centaur");
			case C::Tno:
				return OrbitModel::tr("TNO");
		}
		return {};
	}

	QString periodLabel(orbital_core::CometPeriodClass value)
	{
		using C = orbital_core::CometPeriodClass;
		switch (value) {
			case C::Unknown:
				return {};
			case C::JupiterFamily:
				return OrbitModel::tr("Jupiter family (period)");
			case C::HalleyType:
				return OrbitModel::tr("Halley type (period)");
		}
		return {};
	}

	QString tisserandLabel(orbital_core::CometTisserandClass value)
	{
		using C = orbital_core::CometTisserandClass;
		switch (value) {
			case C::Unknown:
				return {};
			case C::EnckeType:
				return OrbitModel::tr("Encke type (Tj)");
			case C::ChironType:
				return OrbitModel::tr("Chiron type (Tj)");
			case C::JupiterFamily:
				return OrbitModel::tr("Jupiter family (Tj)");
		}
		return {};
	}
}  // namespace

OrbitModel::OrbitModel(QObject* parent) : QObject(parent) {}

OrbitModel::~OrbitModel() = default;

void OrbitModel::setFromFit(const Orbit_p& orbit, const orbital_core::uncertainty::DerivedSigmas& derivedSigmas,
                            bool hasDerivedSigmas, orbital_core::OrbitalClass orbitalClass,
                            orbital_core::CometOrbitClasses cometClasses)
{
	if (!orbit) {
		clear();
		return;
	}

	using namespace Convertors;

	m_orbitType = QString::fromStdString(orbit->getTypeString());
	m_orbitalClass = classLabel(orbitalClass);
	m_orbitalClasses.clear();
	m_orbitalClasses.append(m_orbitalClass.isEmpty() ? m_orbitType : m_orbitalClass);
	for (const auto& label : {periodLabel(cometClasses.period), tisserandLabel(cometClasses.tisserand)}) {
		if (!label.isEmpty())
			m_orbitalClasses.append(label);
	}
	m_e = orbit->getExcentricity();
	m_a = orbit->getSemiMajorAxis();
	m_q = orbit->getPerigee();
	m_i = radiansToDegrees(orbit->getInclination());
	m_node = radiansToDegrees(orbit->getLongitudeOfAscendingNode());
	m_w = radiansToDegrees(orbit->getArgumentOfPericentre());
	m_nu = radiansToDegrees(orbit->getTrueAnomaly());
	m_M = radiansToDegrees(orbit->getMeanAnomaly());
	m_epoch = orbit->getEpoch();

	m_Q = orbit->getApogee();
	m_n = radiansToDegrees(orbit->getMeanMotion());
	m_period = orbit->getOrbitalPeriod();
	m_invA = orbit->getReciprocalSemiMajorAxis();

	m_sigmaE = orbit->getSigmaExcentricity();
	m_sigmaA = orbit->getSigmaSemiMajorAxis();
	m_sigmaI = radiansToDegrees(orbit->getSigmaInclination());
	m_sigmaNode = radiansToDegrees(orbit->getSigmaLongitudeOfAscendingNode());
	m_sigmaW = radiansToDegrees(orbit->getSigmaArgumentOfPericentre());
	m_sigmaNu = radiansToDegrees(orbit->getSigmaTrueAnomaly());
	m_sigmaQ = orbit->getSigmaPerigee();

	m_sigmaM = hasDerivedSigmas ? radiansToDegrees(derivedSigmas.sigma_M) : 0.0;
	m_sigmaPeriod = hasDerivedSigmas ? derivedSigmas.sigma_P : 0.0;

	emit changed();
}

void OrbitModel::clear()
{
	m_orbitType.clear();
	m_orbitalClass.clear();
	m_orbitalClasses.clear();
	m_e = m_a = m_q = m_Q = 0.0;
	m_i = m_node = m_w = m_nu = m_M = 0.0;
	m_n = m_period = m_invA = m_epoch = 0.0;
	m_sigmaE = m_sigmaA = m_sigmaI = m_sigmaNode = 0.0;
	m_sigmaW = m_sigmaNu = m_sigmaM = m_sigmaQ = 0.0;
	m_sigmaPeriod = 0.0;
	emit changed();
}
