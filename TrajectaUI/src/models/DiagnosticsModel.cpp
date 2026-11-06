#include "DiagnosticsModel.h"

#include <QPointF>
#include <QDate>
#include <QXYSeries>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

DiagnosticsModel::DiagnosticsModel(QObject* parent) : QObject(parent) {}

DiagnosticsModel::~DiagnosticsModel() = default;

void DiagnosticsModel::appendSample(const Residual& r, double jd, bool included)
{
	if (r.radar || !std::isfinite(jd) || !std::isfinite(r.dRA_arcsec) || !std::isfinite(r.dDEC_arcsec))
		return;
	m_samples.push_back({jd - 2400000.5, r.dRA_arcsec, r.dDEC_arcsec, included});
}

void DiagnosticsModel::setFromFit(const Residuals& residuals, double optimizationTime, double arcLength)
{
	m_samples.clear();
	m_timeScale = QStringLiteral("TDB");
	for (const auto& r : residuals.data())
		appendSample(r, r.jd, true);
	m_optimizationTime = optimizationTime;
	m_arcLength = arcLength;
	m_rms = residuals.rms();
	emit changed();
}

void DiagnosticsModel::setFromResult(const orbital_core::FitPipelineResult& result, double arcLength)
{
	m_samples.clear();
	m_timeScale = QStringLiteral("UTC");
	const auto& observations = result.observations.all;
	const auto& kept = result.fit.keptObsIndices;
	const auto& all = result.fit.allResiduals.data();
	if (all.size() == observations.size()) {
		std::vector<bool> included(observations.size(), kept.empty());
		for (const auto i : kept)
			if (i < included.size())
				included[i] = true;
		for (std::size_t i = 0; i < all.size(); ++i)
			if (!observations[i].radar)
				appendSample(all[i], observations[i].JD, included[i]);
	} else {
		const auto& fitted = result.fit.residuals.data();
		const auto count = kept.empty() ? observations.size() : kept.size();
		if (fitted.size() == count) {
			for (std::size_t j = 0; j < count; ++j) {
				const auto i = kept.empty() ? j : kept[j];
				if (i < observations.size() && !observations[i].radar)
					appendSample(fitted[j], observations[i].JD, true);
			}
		}
	}
	m_optimizationTime = 0;
	m_arcLength = arcLength;
	m_rms = result.fit.residuals.rms();
	emit changed();
}

void DiagnosticsModel::clear() { setFromFit({}, 0, 0); }

QVariantMap DiagnosticsModel::fillResidualCharts(QObject* draObj, QObject* ddecObj, QObject* excludedDraObj,
                                                 QObject* excludedDdecObj) const
{
	QList<QPointF> ra, dec, excludedRa, excludedDec;
	double start = std::numeric_limits<double>::infinity(), end = -start;
	double raMax = 0.1, decMax = 0.1;
	constexpr int binCount = 12;
	std::array<int, binCount> bins{}, excludedBins{};
	double maxMagnitude = 0;
	for (const auto& sample : m_samples) {
		start = std::min(start, sample.mjd);
		end = std::max(end, sample.mjd);
		maxMagnitude = std::max(maxMagnitude, std::hypot(sample.ra, sample.dec));
	}
	const double origin = m_samples.empty() ? 0 : std::floor(start);
	const double histogramMax = std::max(0.1, std::ceil(maxMagnitude * 10) / 10);
	for (const auto& sample : m_samples) {
		(sample.included ? ra : excludedRa).append({sample.mjd - origin, sample.ra});
		(sample.included ? dec : excludedDec).append({sample.mjd - origin, sample.dec});
		raMax = std::max(raMax, std::abs(sample.ra));
		decMax = std::max(decMax, std::abs(sample.dec));
		const double total = std::hypot(sample.ra, sample.dec);
		const auto bin = std::min(binCount - 1, static_cast<int>(total / histogramMax * binCount));
		++(sample.included ? bins : excludedBins)[bin];
	}
	const auto replace = [](QObject* object, const QList<QPointF>& points) {
		if (auto* series = qobject_cast<QXYSeries*>(object))
			series->replace(points);
	};
	replace(draObj, ra);
	replace(ddecObj, dec);
	replace(excludedDraObj, excludedRa);
	replace(excludedDdecObj, excludedDec);
	QVariantList includedCounts, excludedCounts;
	for (int i = 0; i < binCount; ++i) {
		includedCounts.append(bins[i]);
		excludedCounts.append(excludedBins[i]);
	}
	return {{"count", static_cast<int>(m_samples.size())},
	        {"includedCount", ra.size()},
	        {"excludedCount", excludedRa.size()},
	        {"jdMin", m_samples.empty() ? 0 : start},
	        {"jdMax", m_samples.empty() ? 1 : end},
	        {"draMax", raMax},
	        {"ddecMax", decMax},
	        {"timeOriginMjd", origin},
	        {"timeOriginLabel",
	         m_samples.empty() ? QString()
	                           : QDate::fromJulianDay(static_cast<qint64>(origin) + 2400001).toString(Qt::ISODate) +
	                                 " " + m_timeScale},
	        {"histogramMax", histogramMax},
	        {"bins", includedCounts},
	        {"excludedBins", excludedCounts}};
}
