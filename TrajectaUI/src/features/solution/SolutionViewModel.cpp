#include "SolutionViewModel.h"
#include "util/CalculationTime.h"

#include <orbital_core/CometaryCoordinates.h>
#include <orbital_core/TimeScale.h>
#include <orbital_core/hazard/EquinoctialCoordinates.h>
#include <orbital_core/hazard/ImpactEnergy.h>
#include <orbital_core/hazard/PhaClassifier.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <utility>
#include <QSet>
#include "util/FormatUtils.h"

namespace {
	QString epochText(double jd) { return QString::number(jd, 'f', 9) + QStringLiteral(" JD TDB"); }

	QVariantList equinoctialRows(const orbital_core::FitResult& result)
	{
		namespace equ = orbital_core::hazard::equinoctial;
		namespace cot = orbital_core::cometary;
		const auto& orbit = *result.orbit;
		const auto r = orbit.getR(), v = orbit.getV() * Constants::GAUSS_K;
		equ::Vector6 state;
		state << r.x, r.y, r.z, v.x, v.y, v.z;
		const auto elements = equ::fromCartesian(state);
		if (!elements)
			return {};
		const auto& uncertainty = result.uncertainty;
		std::optional<equ::Matrix6> covariance;
		const auto fromState = equ::cartesianToElementsJacobian(*elements);
		if (fromState) {
			if (result.fitEpochOrbit && orbit.epochMjd() == result.fitEpochOrbit->epochMjd() &&
			    uncertainty.hasStateCovariance()) {
				covariance = (*fromState * uncertainty.cartesianStateCovariance * fromState->transpose()).eval();
			} else {
				cot::Vector6 conventional;
				conventional << orbit.getPerigee(), orbit.getExcentricity(), orbit.getInclination(),
				    orbit.getLongitudeOfAscendingNode(), orbit.getArgumentOfPericentre(), orbit.getTrueAnomaly();
				cot::Matrix6 toState;
				if (cot::toCartesian(conventional, &toState)) {
					const equ::Matrix6 fromCometary = *fromState * toState;
					const auto& cometary = uncertainty.cometaryElementCovariance;
					if (cometary && cometary->epochJdTdb == orbit.getEpoch() && cometary->covariance.allFinite()) {
						covariance = (fromCometary * cometary->covariance * fromCometary.transpose()).eval();
					} else if (uncertainty.hasElementCovariance && uncertainty.elementCovariance.allFinite()) {
						equ::Matrix6 toCometary = equ::Matrix6::Identity();
						toCometary(0, 0) = 1.0 - orbit.getExcentricity();
						toCometary(0, 1) = -orbit.getSemiMajorAxis();
						const equ::Matrix6 fromConventional = fromCometary * toCometary;
						covariance =
						    (fromConventional * uncertainty.elementCovariance * fromConventional.transpose()).eval();
					}
				}
			}
		}
		const QStringList labels{"a", "h", "k", "p", "q", QString::fromUtf8("\xce\xbb")};
		QVariantList rows;
		for (int i = 0; i < 6; ++i) {
			const double scale = i == 5 ? Constants::RAD_TO_DEG : 1.0;
			const double variance = covariance ? (*covariance)(i, i) : std::numeric_limits<double>::quiet_NaN();
			rows.append(QVariantMap{
			    {"label", labels[i]},
			    {"unit", i == 0   ? QStringLiteral("AU")
			             : i == 5 ? QString::fromUtf8("\xc2\xb0")
			                      : QString()},
			    {"value", (*elements)(i)*scale},
			    {"sigma", std::isfinite(variance) && variance >= 0 ? std::sqrt(variance) * scale
			                                                       : std::numeric_limits<double>::quiet_NaN()}});
		}
		return rows;
	}

	QString formatDiameter(double km, bool meters)
	{
		const double value = meters ? km * 1000.0 : km;
		const int digits = meters           ? (value >= 1.0 ? 0 : 2)
		                   : value >= 100.0 ? 0
		                   : value >= 10.0  ? 1
		                   : value >= 1.0   ? 2
		                                    : 3;
		return QString::number(value, 'f', digits);
	}

	QVariantMap sizeRow(double absoluteMagnitude, double albedoLo, double albedoHi)
	{
		const double smallKm = orbital_core::hazard::diameterFromMagnitudeKm(absoluteMagnitude, albedoHi);
		const double largeKm = orbital_core::hazard::diameterFromMagnitudeKm(absoluteMagnitude, albedoLo);
		const bool meters = largeKm < 1.0;
		return {{"label", QString::number(albedoLo, 'f', 2) + QChar(0x2013) + QString::number(albedoHi, 'f', 2)},
		        {"value", formatDiameter(smallKm, meters) + QChar(0x2013) + formatDiameter(largeKm, meters)},
		        {"unit", meters ? QStringLiteral("m") : QStringLiteral("km")},
		        {"albedoLo", albedoLo},
		        {"albedoHi", albedoHi},
		        {"diameterLoKm", smallKm},
		        {"diameterHiKm", largeKm}};
	}
}  // namespace

SolutionViewModel::SolutionViewModel(QObject* parent) : QObject(parent) {}

QVariantMap SolutionViewModel::nominalEarthEntry() const
{
	if (!m_snapshot || !m_snapshot->nominalEarthEntry)
		return {};
	const auto& entry = *m_snapshot->nominalEarthEntry;
	const auto& crossing = entry.crossing;
	const double utc = orbital_core::time::fromMjdTdb(crossing.mjdTdb, orbital_core::time::Scale::UTC);
	const auto coordinate = [](double value, const char* positive, const char* negative) {
		return QString::number(std::abs(value), 'f', 5) + QString::fromUtf8("° ") +
		       QString::fromLatin1(value < 0 ? negative : positive);
	};
	return {{"utc", CalculationTime::formatUtc(utc)},
	        {"latitude", coordinate(crossing.latitudeDeg, "N", "S")},
	        {"longitude", coordinate(crossing.longitudeDeg, "E", "W")},
	        {"altitudeKm", entry.altitudeKm}};
}

QString SolutionViewModel::nominalEarthEntryError() const
{
	return m_snapshot ? QString::fromStdString(m_snapshot->nominalEarthEntryError) : QString();
}

bool SolutionViewModel::publish(const orbital_core::FitPipelineResult& result, const QString& directory,
                                const orbital_core::hazard::PlanetMoidTable& moids)
{
	if (result.fit.cancelled || !result.fit.orbit)
		return false;
	m_snapshot = result;
	m_photometry = orbital_core::photometry::fit(orbital_core::photometry::observations(result.fit).data(),
	                                             result.effectiveConfig.photometry.asteroidG);
	m_cometPhotometry = orbital_core::photometry::fitComet(orbital_core::photometry::observations(result.fit).data(),
	                                                       result.effectiveConfig.photometry);
	m_moids = moids;
	m_equinoctialElements = equinoctialRows(result.fit);
	m_directory = directory;
	const auto& uncertainty = result.fit.uncertainty;
	m_orbit.setFromFit(result.fit.orbit, uncertainty.derivedSigmas,
	                   uncertainty.hasElementCovariance || uncertainty.cometaryElementCovariance.has_value(),
	                   result.orbitalClass(), result.cometClasses());
	m_observations.setResult(result);
	m_covariance.setResult(result.fit);
	m_scene.setOrbit(result.fit.orbit);
	const auto observations = result.fittedObservations();
	m_arcDays = 0;
	if (observations.size() > 1) {
		const auto [lo, hi] = std::ranges::minmax_element(observations, {}, &ObsNode::JD);
		m_arcDays = (*hi).JD - (*lo).JD;
	}
	m_diagnostics.setFromResult(result, m_arcDays);
	emit changed();
	return true;
}

QString SolutionViewModel::designation() const
{
	return m_snapshot ? QString::fromStdString(m_snapshot->observations.objectName) : QString();
}

QString SolutionViewModel::objectTitle() const
{
	if (!m_snapshot)
		return {};
	auto name = QString::fromStdString(Helpers::unpackMpcPackedDesignation(m_snapshot->observations.objectName));
	if (name != designation())
		name.insert(4, ' ');
	return name;
}

QString SolutionViewModel::epoch() const
{
	return m_snapshot ? epochText(m_snapshot->fit.orbit->getEpoch()) : QString();
}

QString SolutionViewModel::status() const
{
	if (!m_snapshot)
		return tr("No solution");
	if (m_snapshot->fit.converged)
		return tr("Converged");
	using Status = orbital_core::FitSelectionStatus;
	switch (m_snapshot->fit.selectionStatus) {
		case Status::PassLimit:
			return tr("Observation selection limit reached");
		case Status::InsufficientObservations:
			return tr("Too few observations after rejection");
		case Status::SafeguardRollback:
			return tr("Observation selection safeguard reached");
		default:
			return tr("Fit not converged");
	}
}

int SolutionViewModel::total() const { return m_snapshot ? static_cast<int>(m_snapshot->observations.all.size()) : 0; }

int SolutionViewModel::included() const
{
	return m_snapshot ? static_cast<int>(m_snapshot->fit.keptObsIndices.size()) : 0;
}

double SolutionViewModel::rms() const { return m_snapshot ? m_snapshot->fit.cost : 0; }

double SolutionViewModel::normalizedRms() const
{
	return m_snapshot ? m_snapshot->fit.normalizedRms : std::numeric_limits<double>::quiet_NaN();
}

bool SolutionViewModel::cometary() const
{
	return orbital_core::cometary::usesCometaryElements(m_orbit.eccentricity());
}

bool SolutionViewModel::hasSigmas() const
{
	return m_snapshot && (m_snapshot->fit.uncertainty.hasElementCovariance ||
	                      m_snapshot->fit.uncertainty.cometaryElementCovariance.has_value());
}

QString SolutionViewModel::epochMjd() const
{
	return m_snapshot ? QString::number(m_snapshot->fit.orbit->epochMjd(), 'f', 9) : QString();
}

QString SolutionViewModel::fitEpochUtc() const
{
	if (!m_snapshot || !m_snapshot->fit.fitEpochOrbit)
		return {};
	const double utc =
	    orbital_core::time::fromMjdTdb(m_snapshot->fit.fitEpochOrbit->epochMjd(), orbital_core::time::Scale::UTC);
	return CalculationTime::formatUtc(utc);
}

QString SolutionViewModel::calculationStartUtc() const
{
	if (!m_snapshot)
		return {};
	const auto utc = CalculationTime::defaultStartUtc(m_snapshot->fit);
	if (!utc)
		return {};
	return CalculationTime::formatUtc(*utc);
}

QString SolutionViewModel::arcRange() const
{
	if (!m_snapshot || !included())
		return {};
	const auto observations = m_snapshot->fittedObservations();
	const auto [lo, hi] = std::ranges::minmax_element(observations, {}, &ObsNode::JD);
	return QString::fromStdString(Convertors::convertJDToGregorian((*lo).JD)).left(10) + " / " +
	       QString::fromStdString(Convertors::convertJDToGregorian((*hi).JD)).left(10);
}

int SolutionViewModel::stationCount() const
{
	QSet<QString> stations;
	if (m_snapshot)
		for (const auto& o : m_snapshot->fittedObservations())
			stations.insert(QString::fromStdString(o.obsCode));
	return static_cast<int>(stations.size());
}

QVariantList SolutionViewModel::stateVector() const
{
	if (!m_snapshot)
		return {};
	const auto r = m_snapshot->fit.orbit->getR();
	const auto v = m_snapshot->fit.orbit->getV() * Constants::GAUSS_K;
	return {r.x, r.y, r.z, v.x, v.y, v.z};
}

bool SolutionViewModel::hasGaiaObservations() const
{
	return m_snapshot &&
	       std::ranges::any_of(m_snapshot->observations.all, [](const auto& obs) { return obs.obsCode == "258"; });
}

bool SolutionViewModel::hasGaiaFprObservations() const
{
	return m_snapshot && std::ranges::any_of(m_snapshot->observations.all, [](const auto& obs) {
		       return obs.astrometricConvention == orbital_core::obs::AstrometricConvention::GaiaFpr;
	       });
}

QString SolutionViewModel::gaiaSummary() const
{
	if (!hasGaiaObservations())
		return {};
	const auto& observations = m_snapshot->observations;
	const auto fprCount = std::ranges::count_if(observations.all, [](const auto& obs) {
		return obs.astrometricConvention == orbital_core::obs::AstrometricConvention::GaiaFpr;
	});
	const auto otherCount = std::ranges::count_if(observations.all, [](const auto& obs) {
		return obs.obsCode == "258" && obs.astrometricConvention != orbital_core::obs::AstrometricConvention::GaiaFpr;
	});
	if (!fprCount)
		return tr("Gaia observations: %1").arg(otherCount);
	QString source;
	switch (observations.gaiaRadiusSource) {
		case orbital_core::GaiaRadiusSource::Sbdb:
			source = tr("SBDB");
			break;
		case orbital_core::GaiaRadiusSource::FitEstimate:
			source = tr("fit estimate");
			break;
		case orbital_core::GaiaRadiusSource::None:
			source = tr("no correction");
			break;
		case orbital_core::GaiaRadiusSource::Unknown:
			break;
	}
	const auto radius = observations.gaiaPhotocenterRadiusKm.value_or(0.0);
	const auto summary = source.isEmpty()
	                         ? tr("Gaia FPR: %1 · radius %2 km").arg(fprCount).arg(radius, 0, 'g', 5)
	                         : tr("Gaia FPR: %1 · radius %2 km (%3)").arg(fprCount).arg(radius, 0, 'g', 5).arg(source);
	return otherCount ? tr("Gaia observations: %1 · %2").arg(otherCount).arg(summary) : summary;
}

double SolutionViewModel::absoluteMagnitude() const
{
	if (cometPhotometry())
		return std::numeric_limits<double>::quiet_NaN();
	return m_photometry.absoluteMagnitude.value_or(std::numeric_limits<double>::quiet_NaN());
}

bool SolutionViewModel::cometPhotometry() const
{
	if (!m_snapshot)
		return false;
	return m_snapshot->usesCometPhotometry();
}

double SolutionViewModel::cometTotalMagnitude() const
{
	return m_cometPhotometry.total.absoluteMagnitude.value_or(std::numeric_limits<double>::quiet_NaN());
}

double SolutionViewModel::cometNuclearMagnitude() const
{
	return m_cometPhotometry.nuclear.absoluteMagnitude.value_or(std::numeric_limits<double>::quiet_NaN());
}

QVariantList SolutionViewModel::sizeEstimates() const
{
	const double h = absoluteMagnitude();
	if (!m_snapshot || !std::isfinite(h))
		return {};
	constexpr std::array<std::pair<double, double>, 5> bands{
	    {{0.05, 0.10}, {0.10, 0.15}, {0.15, 0.25}, {0.25, 0.50}, {0.50, 1.00}}};
	QVariantList rows;
	for (const auto& [lo, hi] : bands)
		rows.append(sizeRow(h, lo, hi));
	return rows;
}

QVariantMap SolutionViewModel::typicalSize() const
{
	const double h = absoluteMagnitude();
	return m_snapshot && std::isfinite(h) ? sizeRow(h, 0.05, 0.25) : QVariantMap{};
}

double SolutionViewModel::earthMoidKm() const
{
	return m_moids[orbital_core::hazard::Planet::Earth] * (Constants::AU_IN_METERS / 1000.0);
}

orbital_core::hazard::PhaStatus SolutionViewModel::phaClassification() const
{
	if (!m_snapshot)
		return orbital_core::hazard::PhaStatus::NotApplicable;
	using namespace orbital_core::hazard;
	return classifyPha(!cometPhotometry(), absoluteMagnitude(), m_moids[Planet::Earth]);
}

bool SolutionViewModel::potentiallyHazardous() const
{
	return phaClassification() == orbital_core::hazard::PhaStatus::Pha;
}

QString SolutionViewModel::phaStatus() const
{
	using orbital_core::hazard::PhaStatus;
	switch (phaClassification()) {
		case PhaStatus::Pha:
			return tr("PHA");
		case PhaStatus::NotPha:
			return tr("Not PHA");
		case PhaStatus::Unknown:
			return tr("Unknown");
		case PhaStatus::NotApplicable:
			return {};
	}
	return {};
}

QVariantList SolutionViewModel::planetMoids() const
{
	QVariantList rows;
	for (const auto planet : orbital_core::hazard::kAllPlanets) {
		if (planet == orbital_core::hazard::Planet::Moon)
			continue;
		const auto label = orbital_core::hazard::longLabel(planet);
		const double au = m_moids[planet];
		rows.append(QVariantMap{{"planet", QString::fromUtf8(label.data(), label.size())},
		                        {"au", au},
		                        {"km", au * (Constants::AU_IN_METERS / 1000.0)}});
	}
	return rows;
}

bool SolutionViewModel::nonGravEnabled() const { return m_snapshot && m_snapshot->fit.finalNonGrav.enabled(); }

QString SolutionViewModel::nonGravModel() const
{
	if (!m_snapshot)
		return {};
	if (!nonGravEnabled())
		return tr("Off");
	switch (m_snapshot->fit.finalNonGrav.model()) {
		case NgModel::MarsdenH2O:
			return tr("Marsden H2O");
		case NgModel::MarsdenCO:
			return tr("Marsden CO");
		case NgModel::InverseSquare:
			return tr("Inverse-square");
		case NgModel::Yarkovsky:
			return tr("Yarkovsky");
			return tr("Inverse-square");
		case NgModel::None:
			return tr("Off");
	}
	return tr("Unknown");
}

QVariantList SolutionViewModel::nonGravParameters() const
{
	if (!m_snapshot)
		return {};
	const auto& result = m_snapshot->fit;
	const auto& model = result.finalNonGrav;
	const auto supported =
	    model.enabled() ? NonGravModel::defaultFitParameterMaskFor(model.model()) : NgParameterMask::None;
	const double missing = std::numeric_limits<double>::quiet_NaN();
	const std::array values{model.A1(), model.A2(), model.A3(), model.DT()};
	std::array sigmas{missing, missing, missing, missing};
	const auto& covariance = result.uncertainty.fitParameterCovariance;
	const int dimension = 6 + model.numFitParams();
	if (result.fitEpochOrbit && covariance.rows() == dimension && covariance.cols() == dimension &&
	    covariance.allFinite() && !covariance.isZero(0)) {
		for (int i = 0; i < model.numFitParams(); ++i) {
			const double variance = covariance(6 + i, 6 + i);
			if (variance >= 0)
				sigmas[static_cast<int>(model.fitParameterAt(i))] = std::sqrt(variance);
		}
	}
	const QStringList labels{"A1", "A2", "A3", "DT"};
	QVariantList rows;
	for (int i = 0; i < NonGravModel::maxFitParams; ++i) {
		const auto bit = static_cast<NgParameterMask>(1u << i);
		const bool estimated = hasNgParameter(model.fitParameterMask(), bit);
		const bool used = model.enabled() && (hasNgParameter(supported, bit) || estimated ||
		                                      (values[i] != 0.0 && (i < 3 || model.model() == NgModel::InverseSquare)));
		rows.append(QVariantMap{{"label", labels[i]},
		                        {"unit", i == 3 ? tr("days") : QStringLiteral("AU/day\u00b2")},
		                        {"value", used ? values[i] : missing},
		                        {"sigma", sigmas[i]},
		                        {"used", used},
		                        {"estimated", estimated},
		                        {"status", !used       ? tr("Not used")
		                                   : estimated ? tr("Estimated")
		                                               : tr("Fixed")}});
	}
	return rows;
}
