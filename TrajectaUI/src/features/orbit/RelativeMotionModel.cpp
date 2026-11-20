#include "RelativeMotionModel.h"
#include "RelativePathSegments.h"
#include <orbital_core/hazard/PlanetData.h>
#include <orbital_core/FitPipeline.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <QDate>
#include <QMetaObject>
#include "util/FormatUtils.h"
#include "util/CalculationTime.h"

RelativeMotionModel::RelativeMotionModel(QObject* parent) : QObject(parent) {}

RelativeMotionModel::~RelativeMotionModel()
{
	m_relativeWorker.request_stop();
	if (m_relativeWorker.joinable())
		m_relativeWorker.join();
}

void RelativeMotionModel::reset()
{
	m_relativeWorker.request_stop();
	++m_relativeGeneration;
	m_relativeBusy = false;
	m_relativeSolution.reset();
	m_relativeSamples.clear();
	m_captureIntervals.clear();
	m_relativePath.clear();
	m_relativeMoonPath.clear();
	m_relativeBoundSegments.clear();
	m_relativeUnboundSegments.clear();
	m_relativeDistantNegativeSegments.clear();
	m_relativeCaptures.clear();
	m_relativeCaptureMarkers.clear();
	m_relativeError.clear();
	m_position = 0;
	m_objectPosition = {};
	m_relativeMoon = {};
	emit relativeChanged();
	emit statusChanged();
}

void RelativeMotionModel::setPosition(double position)
{
	m_position = position;
	updateRelativePosition();
	emit positionChanged();
}

std::optional<double> RelativeMotionModel::captureStartEpoch(int index) const
{
	if (index < 0 || static_cast<std::size_t>(index) >= m_captureIntervals.size())
		return std::nullopt;
	return m_captureIntervals[static_cast<std::size_t>(index)].startMjdTdb;
}

void RelativeMotionModel::setSolution(const orbital_core::FitPipelineResult& result)
{
	const auto& orbit = result.fit.fitEpochOrbit ? result.fit.fitEpochOrbit : result.fit.orbit;
	if (!orbit || !result.fit.converged)
		return;
	orbital_core::EphemerisSolution solution;
	solution.position = orbit->getR();
	solution.velocity = orbit->getV() * Constants::GAUSS_K;
	solution.epochMjdTdb = orbit->epochMjd();
	solution.forces = result.effectiveConfig.nbody;
	solution.nonGrav = result.fit.finalNonGrav;
	if (result.fit.fitEpochOrbit)
		solution.covariance = result.fit.uncertainty.fitParameterCovariance;
	m_relativeSolution = solution;
	const double utc = orbital_core::time::fromMjdTdb(solution.epochMjdTdb, orbital_core::time::Scale::UTC);
	if (std::isfinite(utc)) {
		const QDate date = QDate::fromJulianDay(static_cast<qint64>(std::floor(utc)) + 2400001);
		m_relativeStartDate = date.addDays(-180).toString(Qt::ISODate);
		m_relativeEndDate = date.addDays(365).toString(Qt::ISODate);
	}
	emit settingsChanged();
}

void RelativeMotionModel::setRelativeTimeScale(const QString& value)
{
	if ((value != QLatin1String("UTC") && value != QLatin1String("TDB")) || m_relativeTimeScale == value)
		return;
	m_relativeTimeScale = value;
	emit settingsChanged();
}

void RelativeMotionModel::setRotatingFrame(bool value)
{
	if (m_rotatingFrame == value)
		return;
	m_rotatingFrame = value;
	rebuildRelativePath();
	updateRelativePosition();
	emit relativeChanged();
	emit positionChanged();
}

void RelativeMotionModel::setRelativeStartDate(const QString& value)
{
	if (m_relativeStartDate == value)
		return;
	m_relativeStartDate = value;
	emit settingsChanged();
}

void RelativeMotionModel::setRelativeEndDate(const QString& value)
{
	if (m_relativeEndDate == value)
		return;
	m_relativeEndDate = value;
	emit settingsChanged();
}

void RelativeMotionModel::setRelativeStepDays(double value)
{
	if (m_relativeStepDays == value || (std::isnan(m_relativeStepDays) && std::isnan(value)))
		return;
	m_relativeStepDays = value;
	emit settingsChanged();
}

void RelativeMotionModel::setRelativePlanetIndex(int value)
{
	if (m_relativePlanetIndex == value)
		return;
	m_relativePlanetIndex = value;
	emit settingsChanged();
}

void RelativeMotionModel::rebuildRelativePath()
{
	m_relativePath.clear();
	m_relativeMoonPath.clear();
	m_relativeBoundSegments.clear();
	m_relativeUnboundSegments.clear();
	m_relativeDistantNegativeSegments.clear();
	m_relativeRadius = 1e-5;
	for (const auto& sample : m_relativeSamples) {
		const auto r = m_rotatingFrame ? sample.rotatingPositionAu : sample.positionAu;
		m_relativePath.append(QVariant::fromValue(FormatUtils::eclToQt3D(r)));
		if (sample.hasMoon)
			m_relativeMoonPath.append(QVariant::fromValue(
			    FormatUtils::eclToQt3D(m_rotatingFrame ? sample.rotatingMoonPositionAu : sample.moonPositionAu)));
		m_relativeRadius = std::max({m_relativeRadius, r.magnitude(), sample.hillRadiusAu});
	}
	const auto appendSegment = [](QVariantList& target, const QVector3D& a, const QVector3D& b) {
		target.append(QVariant::fromValue(a));
		target.append(QVariant::fromValue(b));
	};
	for (std::size_t i = 1; i < m_relativeSamples.size(); ++i) {
		const auto& a = m_relativeSamples[i - 1];
		const auto& b = m_relativeSamples[i];
		const auto p = m_relativePath.at(static_cast<qsizetype>(i - 1)).value<QVector3D>();
		const auto q = m_relativePath.at(static_cast<qsizetype>(i)).value<QVector3D>();
		const auto parts = relative_path::split(a.specificEnergyAu2PerDay2, b.specificEnergyAu2PerDay2,
		                                        a.distanceAu - 3 * a.hillRadiusAu, b.distanceAu - 3 * b.hillRadiusAu);
		for (std::size_t j = 0; j < parts.count; ++j) {
			const auto& part = parts.values[j];
			auto& target = part.kind == relative_path::Kind::Unbound   ? m_relativeUnboundSegments
			               : part.kind == relative_path::Kind::Capture ? m_relativeBoundSegments
			                                                           : m_relativeDistantNegativeSegments;
			appendSegment(target, p + (q - p) * static_cast<float>(part.start),
			              p + (q - p) * static_cast<float>(part.end));
		}
	}
	m_relativeCaptureMarkers.clear();
	const auto addMarker = [&](double epoch, int eventIndex, bool isStart) {
		const auto next = std::lower_bound(m_relativeSamples.begin(), m_relativeSamples.end(), epoch,
		                                   [](const auto& sample, double t) { return sample.mjdTdb < t; });
		const auto last = next == m_relativeSamples.end() ? std::prev(m_relativeSamples.end()) : next;
		const auto first = last == m_relativeSamples.begin() ? last : std::prev(last);
		const double fraction = first == last ? 0 : (epoch - first->mjdTdb) / (last->mjdTdb - first->mjdTdb);
		const auto a = m_rotatingFrame ? first->rotatingPositionAu : first->positionAu;
		const auto b = m_rotatingFrame ? last->rotatingPositionAu : last->positionAu;
		QVariantMap marker;
		marker.insert("position", QVariant::fromValue(FormatUtils::eclToQt3D(a * (1 - fraction) + b * fraction)));
		marker.insert("event", eventIndex + 1);
		marker.insert("start", isStart);
		m_relativeCaptureMarkers.append(marker);
	};
	if (!m_relativeSamples.empty()) {
		for (std::size_t i = 0; i < m_captureIntervals.size(); ++i) {
			const auto& event = m_captureIntervals[i];
			if (!event.startClipped)
				addMarker(event.startMjdTdb, static_cast<int>(i), true);
			if (!event.endClipped)
				addMarker(event.endMjdTdb, static_cast<int>(i), false);
		}
	}
}

const orbital_core::RelativeMotionSample* RelativeMotionModel::currentRelativeSample() const
{
	if (!relativeAvailable())
		return nullptr;
	const double epoch = m_relativeSamples.front().mjdTdb +
	                     m_position * (m_relativeSamples.back().mjdTdb - m_relativeSamples.front().mjdTdb);
	const auto next = std::lower_bound(m_relativeSamples.begin(), m_relativeSamples.end(), epoch,
	                                   [](const auto& sample, double time) { return sample.mjdTdb < time; });
	if (next == m_relativeSamples.begin())
		return &*next;
	if (next == m_relativeSamples.end())
		return &m_relativeSamples.back();
	const auto previous = std::prev(next);
	return epoch - previous->mjdTdb <= next->mjdTdb - epoch ? &*previous : &*next;
}

void RelativeMotionModel::updateRelativePosition()
{
	if (!relativeAvailable())
		return;
	const double epoch = m_relativeSamples.front().mjdTdb +
	                     m_position * (m_relativeSamples.back().mjdTdb - m_relativeSamples.front().mjdTdb);
	const auto next = std::lower_bound(m_relativeSamples.begin(), m_relativeSamples.end(), epoch,
	                                   [](const auto& sample, double time) { return sample.mjdTdb < time; });
	const auto last = next == m_relativeSamples.end() ? std::prev(m_relativeSamples.end()) : next;
	const auto first = last == m_relativeSamples.begin() ? last : std::prev(last);
	const double fraction = first == last ? 0 : (epoch - first->mjdTdb) / (last->mjdTdb - first->mjdTdb);
	const auto lerp = [fraction](const Vec3D& a, const Vec3D& b) { return a * (1 - fraction) + b * fraction; };
	const auto& a = *first;
	const auto& b = *last;
	m_objectPosition = FormatUtils::eclToQt3D(lerp(m_rotatingFrame ? a.rotatingPositionAu : a.positionAu,
	                                               m_rotatingFrame ? b.rotatingPositionAu : b.positionAu));
	m_relativeMoon = FormatUtils::eclToQt3D(lerp(m_rotatingFrame ? a.rotatingMoonPositionAu : a.moonPositionAu,
	                                             m_rotatingFrame ? b.rotatingMoonPositionAu : b.moonPositionAu));
}

double RelativeMotionModel::relativeHillRadius() const
{
	const auto* s = currentRelativeSample();
	return s ? s->hillRadiusAu : 0;
}

double RelativeMotionModel::relativeDistanceKm() const
{
	const auto* s = currentRelativeSample();
	return s ? s->distanceAu * Constants::AU_IN_METERS / 1000 : 0;
}

double RelativeMotionModel::relativeEnergy() const
{
	const auto* s = currentRelativeSample();
	return s ? s->specificEnergyAu2PerDay2 * std::pow(Constants::AU_IN_METERS / 1000 / Constants::SECONDS_PER_DAY, 2)
	         : 0;
}

QVector3D RelativeMotionModel::relativeSunDirection() const
{
	const auto* s = currentRelativeSample();
	return s ? FormatUtils::eclToQt3D(m_rotatingFrame ? Vec3D(-1, 0, 0) : s->sunDirection) : QVector3D();
}

void RelativeMotionModel::calculateRelative(const QString& startUtc, const QString& endUtc, int planetIndex,
                                            double stepDays)
{
	if (m_relativeBusy)
		return;
	const QDate start = QDate::fromString(startUtc.trimmed(), Qt::ISODate);
	const QDate end = QDate::fromString(endUtc.trimmed(), Qt::ISODate);
	if (!m_relativeSolution || !start.isValid() || !end.isValid() || end <= start || planetIndex < 0 ||
	    planetIndex > 7 || !std::isfinite(stepDays) || stepDays <= 0) {
		m_relativeError = tr("Choose a valid interval, planet and step.");
		emit statusChanged();
		return;
	}
	if (m_relativeWorker.joinable())
		m_relativeWorker.join();
	orbital_core::RelativeMotionRequest request;
	request.solution = *m_relativeSolution;
	request.planet = static_cast<EphemerisPlanet>(static_cast<int>(EphemerisPlanet::Mercury) + planetIndex);
	const double first = start.toJulianDay() - 2400001.0;
	const double last = end.toJulianDay() - 2400001.0;
	const QString selectedTimeScale = m_relativeTimeScale;
	const auto scale =
	    selectedTimeScale == QLatin1String("TDB") ? orbital_core::time::Scale::TDB : orbital_core::time::Scale::UTC;
	if (!std::isfinite(orbital_core::time::toMjdTdb(first, scale)) ||
	    !std::isfinite(orbital_core::time::toMjdTdb(last, scale))) {
		m_relativeError = tr("UTC is unavailable before 1960. Select TDB for historical dates.");
		emit statusChanged();
		return;
	}
	constexpr std::size_t maxPoints = 100000;
	const double estimated = std::ceil((last - first) / stepDays) + 1;
	if (!std::isfinite(estimated) || first + stepDays <= first) {
		m_relativeError = tr("The step is too small to represent this trajectory.");
		emit statusChanged();
		return;
	}
	if (estimated > static_cast<double>(maxPoints)) {
		m_relativeError = tr("The trajectory exceeds the %1-point limit. Increase the step or shorten the interval.")
		                      .arg(maxPoints);
		emit statusChanged();
		return;
	}
	request.maxRefinedSamples = maxPoints;
	const quint64 generation = ++m_relativeGeneration;
	m_relativeBusy = true;
	m_relativeError.clear();
	emit statusChanged();
	const QString center =
	    QStringList{"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"}.at(planetIndex);
	m_relativeWorker = std::jthread([this, generation, center, planetIndex, first, last, stepDays, estimated, scale,
	                                 selectedTimeScale, request = std::move(request)](std::stop_token stop) mutable {
		try {
			request.epochsMjdTdb.reserve(static_cast<std::size_t>(estimated));
			for (std::size_t i = 0; i < static_cast<std::size_t>(estimated) && !stop.stop_requested(); ++i) {
				const double epoch = first + i * stepDays;
				if (epoch >= last)
					break;
				request.epochsMjdTdb.push_back(orbital_core::time::toMjdTdb(epoch, scale));
			}
			auto samples = std::make_shared<std::vector<orbital_core::RelativeMotionSample>>();
			auto captures = std::make_shared<std::vector<orbital_core::TemporaryCaptureInterval>>();
			auto hillPassages = std::make_shared<std::vector<orbital_core::HillSpherePassage>>();
			auto uncertainties = std::make_shared<orbital_core::RelativeMotionTimeUncertainty>();
			if (!stop.stop_requested()) {
				request.epochsMjdTdb.push_back(orbital_core::time::toMjdTdb(last, scale));
				*samples = orbital_core::computeRefinedRelativeMotion(request, {.cancel = stop});
				if (!stop.stop_requested()) {
					*captures = orbital_core::findTemporaryCaptures(*samples);
					*hillPassages = orbital_core::findHillSpherePassages(*samples);
					try {
						*uncertainties = orbital_core::relativeMotionTimeUncertainty(request, *samples, *captures,
						                                                             *hillPassages, {.cancel = stop});
					} catch (const std::exception&) {
						*uncertainties = {};
					}
				}
			}
			QMetaObject::invokeMethod(
			    this,
			    [this, generation, center, planetIndex, samples, captures, hillPassages, uncertainties, scale,
			     selectedTimeScale] {
				    if (generation != m_relativeGeneration)
					    return;
				    m_relativeBusy = false;
				    const bool accepted = !samples->empty();
				    if (samples->empty()) {
					    m_relativeError = tr("Calculation cancelled.");
				    } else {
					    m_relativeSamples = std::move(*samples);
					    m_relativeResultTimeScale = selectedTimeScale;
					    m_relativeCenter = center;
					    m_relativePlanetRadius = orbital_core::hazard::planetRadiusKm(
					                                 static_cast<orbital_core::hazard::Planet>(planetIndex)) /
					                             (Constants::AU_IN_METERS / 1000.0);
					    m_captureIntervals = std::move(*captures);
					    m_relativeCaptures.clear();
					    for (std::size_t i = 0; i < m_captureIntervals.size(); ++i) {
						    const auto& event = m_captureIntervals[i];
						    QVariantMap row;
						    row.insert("start", CalculationTime::formatUtc(
						                            orbital_core::time::fromMjdTdb(event.startMjdTdb, scale)));
						    row.insert("end", CalculationTime::formatUtc(
						                          orbital_core::time::fromMjdTdb(event.endMjdTdb, scale)));
						    row.insert("durationDays", event.endMjdTdb - event.startMjdTdb);
						    row.insert("startMjdTdb", event.startMjdTdb);
						    row.insert("endMjdTdb", event.endMjdTdb);
						    row.insert("revolutions", event.signedRevolutions);
						    row.insert("angularTravel", event.angularTravelRevolutions);
						    row.insert("revolutionsValid", event.revolutionsValid);
						    row.insert("inertialRevolutions", event.inertialSignedRevolutions);
						    row.insert("inertialAngularTravel", event.inertialAngularTravelRevolutions);
						    row.insert("inertialRevolutionsValid", event.inertialRevolutionsValid);
						    row.insert("startClipped", event.startClipped);
						    row.insert("endClipped", event.endClipped);
						    QVariantList passages;
						    for (std::size_t j = 0; j < hillPassages->size(); ++j) {
							    const auto& passage = (*hillPassages)[j];
							    if (passage.exitMjdTdb <= event.startMjdTdb || passage.entryMjdTdb >= event.endMjdTdb)
								    continue;
							    QVariantMap crossing;
							    crossing.insert("entry", CalculationTime::formatUtc(orbital_core::time::fromMjdTdb(
							                                 passage.entryMjdTdb, scale)));
							    crossing.insert("exit", CalculationTime::formatUtc(
							                                orbital_core::time::fromMjdTdb(passage.exitMjdTdb, scale)));
							    crossing.insert("entryClipped", passage.entryClipped);
							    crossing.insert("exitClipped", passage.exitClipped);
							    if (j < uncertainties->hillPassages.size()) {
								    const auto& uncertainty = uncertainties->hillPassages[j];
								    if (uncertainty.start3SigmaDays)
									    crossing.insert("entry3SigmaSeconds", *uncertainty.start3SigmaDays * 86400);
								    if (uncertainty.end3SigmaDays)
									    crossing.insert("exit3SigmaSeconds", *uncertainty.end3SigmaDays * 86400);
							    }
							    passages.append(crossing);
						    }
						    row.insert("hillPassages", passages);
						    if (i < uncertainties->captures.size()) {
							    const auto& uncertainty = uncertainties->captures[i];
							    if (uncertainty.start3SigmaDays)
								    row.insert("start3SigmaSeconds", *uncertainty.start3SigmaDays * 86400);
							    if (uncertainty.end3SigmaDays)
								    row.insert("end3SigmaSeconds", *uncertainty.end3SigmaDays * 86400);
						    }
						    m_relativeCaptures.append(row);
					    }
					    m_position = 0;
					    rebuildRelativePath();
					    updateRelativePosition();
				    }
				    if (accepted)
					    emit relativeChanged();
				    emit statusChanged();
				    emit resultReady(accepted);
				    emit positionChanged();
			    },
			    Qt::QueuedConnection);
		} catch (const std::exception& e) {
			const QString message = QString::fromUtf8(e.what());
			QMetaObject::invokeMethod(
			    this,
			    [this, generation, message] {
				    if (generation != m_relativeGeneration)
					    return;
				    m_relativeBusy = false;
				    m_relativeError = message;
				    emit statusChanged();
			    },
			    Qt::QueuedConnection);
		}
	});
}

void RelativeMotionModel::cancelRelative()
{
	if (m_relativeBusy)
		m_relativeWorker.request_stop();
}
