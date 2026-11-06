#include "ObservationTableModel.h"

#include <QAbstractTableModel>
#include <QSet>
#include <QtCharts/QXYSeries>
#include "util/FormatUtils.h"
#include "orbital_core/obs/StarCatalog.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <limits>
#include <unordered_map>

namespace {
	constexpr int numericRole = Qt::UserRole + 1;
	constexpr int sourceRowRole = Qt::UserRole + 2;
	constexpr int includedRole = Qt::UserRole + 3;
	constexpr int desiredIncludedRole = Qt::UserRole + 4;
	constexpr int pendingRole = Qt::UserRole + 5;
	constexpr int manualRole = Qt::UserRole + 6;
	using Disposition = orbital_core::ObservationDisposition;

	double rangeRateScale(const orbital_core::obs::RadarObservation& radar)
	{
		return -Constants::AU_PER_DAY * Constants::AU_IN_METERS / 86400.0 / (2.0 * radar.frequencyHz);
	}

	QString number(double value, int precision = 3)
	{
		return std::isfinite(value) ? QString::number(value, 'f', precision) : QStringLiteral("-");
	}

	QString declination(double rad) { return FormatUtils::dec(rad); }

	QString timeScaleName(orbital_core::time::Scale scale)
	{
		using orbital_core::time::Scale;
		switch (scale) {
			case Scale::UTC:
				return QStringLiteral("UTC");
			case Scale::UT1:
				return QStringLiteral("UT1");
			case Scale::TT:
				return QStringLiteral("TT");
			case Scale::TDB:
				return QStringLiteral("TDB");
		}
		return {};
	}

	QString techniqueName(char code)
	{
		static const QHash<char, QString> names{{'C', "CCD"},
		                                        {'c', "Corrected CCD"},
		                                        {'B', "CMOS"},
		                                        {'P', "Photographic"},
		                                        {'A', "Precessed from B1950"},
		                                        {'H', "Photometer"},
		                                        {'X', "Replaced observation"},
		                                        {'x', "Replaced observation"},
		                                        {'n', "Video"},
		                                        {'e', "Encoder"},
		                                        {'M', "Micrometer"},
		                                        {'T', "Transit circle"},
		                                        {'E', "Occultation"},
		                                        {'N', "Normal place"}};
		if (code == ' ' || code == '\0')
			return QStringLiteral("Not supplied");
		return names.contains(code) ? names.value(code) : QStringLiteral("MPC code %1").arg(QChar::fromLatin1(code));
	}

	QString catalogName(const ObsNode& observation)
	{
		const char code = observation.catalogCode;
		const bool noCode = code == ' ' || code == '\0';
		if (observation.astCat.empty() && noCode)
			return QStringLiteral("Not supplied");
		const auto* catalog = observation.astCat.empty() ? orbital_core::obs::findStarCatalog(code)
		                                                 : orbital_core::obs::findStarCatalog(observation.astCat);
		if (!catalog)
			return observation.astCat.empty()
			           ? QStringLiteral("MPC code %1").arg(QChar::fromLatin1(code))
			           : QString::fromStdString(observation.astCat) + QStringLiteral(" (unrecognized)");
		const auto name = QString::fromUtf8(catalog->name.data(), static_cast<qsizetype>(catalog->name.size()));
		return noCode ? name : name + " (" + QChar::fromLatin1(code) + ")";
	}
}  // namespace

class ObservationRows final : public QAbstractTableModel {
   public:
	struct Row {
		orbital_core::ObservationId id;
		Disposition accepted{Disposition::Auto};
		Disposition draft{Disposition::Auto};
		QString epoch;
		QString timeScale;
		QString station;
		QString technique;
		QString catalog;
		QString program;
		double jd{};
		double ra{};
		double dec{};
		bool included{};
		bool replaced{};
		bool debiasEnabled{};
		double biasRaArcsec{};
		double biasDecArcsec{};
		std::optional<Residual> residual;
		std::optional<orbital_core::obs::RadarObservation> radar;

		bool desiredIncluded() const { return draft == Disposition::Auto ? included : draft == Disposition::Include; }

		bool pending() const { return draft != accepted; }
	};

	using QAbstractTableModel::QAbstractTableModel;
	std::vector<Row> rows;
	bool radarMode{};
	bool rangeRates{};

	int rowCount(const QModelIndex& parent = {}) const override
	{
		return parent.isValid() ? 0 : static_cast<int>(rows.size());
	}

	int columnCount(const QModelIndex& parent = {}) const override { return parent.isValid() ? 0 : 11; }

	QHash<int, QByteArray> roleNames() const override
	{
		return {{Qt::DisplayRole, "display"},
		        {numericRole, "numeric"},
		        {sourceRowRole, "sourceRow"},
		        {includedRole, "included"},
		        {desiredIncludedRole, "desiredIncluded"},
		        {pendingRole, "pending"},
		        {manualRole, "manual"}};
	}

	QVariant headerData(int section, Qt::Orientation orientation, int role) const override
	{
		if (role != Qt::DisplayRole)
			return {};
		if (orientation == Qt::Vertical)
			return section + 1;
		if (radarMode) {
			const QStringList labels{tr("No."),
			                         tr("Date and time"),
			                         tr("Receiver"),
			                         tr("Delay, s"),
			                         rangeRates ? tr("Range rate, m/s") : tr("Doppler, Hz"),
			                         tr("O-C delay, us"),
			                         rangeRates ? tr("O-C range rate, m/s") : tr("O-C Doppler, Hz"),
			                         tr("\u03c3 delay, us"),
			                         rangeRates ? tr("\u03c3 range rate, m/s") : tr("\u03c3 Doppler, Hz"),
			                         QStringLiteral("\u03c7 ") + tr("(meas.)"),
			                         tr("Fit")};
			return labels.value(section);
		}
		const QStringList labels{tr("No."),
		                         tr("Date and time"),
		                         tr("Station"),
		                         tr("RA"),
		                         tr("Dec"),
		                         QStringLiteral("\u0394\u03b1 cos \u03b4, \u2033"),
		                         QStringLiteral("\u0394\u03b4, \u2033"),
		                         QStringLiteral("\u03c3 RA cos \u03b4, \u2033"),
		                         QStringLiteral("\u03c3 Dec, \u2033"),
		                         QStringLiteral("\u03c7 ") + tr("(meas.)"),
		                         tr("Fit")};
		return labels.value(section);
	}

	QVariant data(const QModelIndex& index, int role) const override
	{
		if (!index.isValid() || index.row() < 0 || index.row() >= rowCount() || index.column() >= columnCount())
			return {};
		const auto& r = rows[index.row()];
		if (role == sourceRowRole)
			return index.row();
		if (role == includedRole)
			return r.included;
		if (role == desiredIncludedRole)
			return r.desiredIncluded();
		if (role == pendingRole)
			return r.pending();
		if (role == manualRole)
			return r.draft != Disposition::Auto;
		if (role != Qt::DisplayRole && role != numericRole)
			return {};
		const bool numeric = role == numericRole;
		if (r.radar && index.column() >= 3 && index.column() <= 9) {
			const auto* residual = r.residual && r.residual->radar ? &*r.residual->radar : nullptr;
			std::optional<double> value;
			int precision = 3;
			switch (index.column()) {
				case 3:
					if (r.radar->delaySeconds)
						value = r.radar->delaySeconds->value;
					precision = 8;
					break;
				case 4:
					if (r.radar->dopplerHz)
						value = r.radar->dopplerHz->value;
					break;
				case 5:
					if (residual && residual->delaySeconds)
						value = residual->delaySeconds->residual * 1e6;
					break;
				case 6:
					if (residual && residual->dopplerHz)
						value = residual->dopplerHz->residual;
					break;
				case 7:
					if (r.radar->delaySeconds)
						value = r.radar->delaySeconds->sigma * 1e6;
					break;
				case 8:
					if (r.radar->dopplerHz)
						value = r.radar->dopplerHz->sigma;
					break;
				case 9:
					if (residual)
						value = std::sqrt(residual->chiSquared);
					break;
			}
			if (value && rangeRates && (index.column() == 4 || index.column() == 6 || index.column() == 8)) {
				const double scale = rangeRateScale(*r.radar);
				*value *= index.column() == 8 ? std::abs(scale) : scale;
				precision = 6;
			}
			return !value ? QVariant() : numeric ? QVariant(*value) : QVariant(number(*value, precision));
		}
		const auto* weight = r.residual && r.residual->astrometricWeight ? &*r.residual->astrometricWeight : nullptr;
		const auto weightedNumber = [&](double value, int precision) -> QVariant {
			if ((r.replaced && !r.included) || !weight || !weight->weighted || !std::isfinite(value) || value < 0)
				return {};
			return numeric ? QVariant(value) : QVariant(number(value, precision));
		};
		switch (index.column()) {
			case 0:
				return index.row() + 1;
			case 1:
				return numeric ? QVariant(r.jd) : QVariant(r.epoch + QLatin1Char(' ') + r.timeScale);
			case 2:
				return r.station;
			case 3:
				return numeric ? QVariant(r.ra) : QVariant(FormatUtils::ra(r.ra));
			case 4:
				return numeric ? QVariant(r.dec) : QVariant(declination(r.dec));
			case 5:
				return r.residual
				           ? (numeric ? QVariant(r.residual->dRA_arcsec) : QVariant(number(r.residual->dRA_arcsec)))
				           : QVariant();
			case 6:
				return r.residual
				           ? (numeric ? QVariant(r.residual->dDEC_arcsec) : QVariant(number(r.residual->dDEC_arcsec)))
				           : QVariant();
			case 7:
				return weightedNumber(weight ? weight->sigmaRaArcsec : -1, 3);
			case 8:
				return weightedNumber(weight ? weight->sigmaDecArcsec : -1, 3);
			case 9:
				return weightedNumber(weight && weight->chiSquared >= 0 ? std::sqrt(weight->chiSquared) : -1, 2);
			case 10:
				if (numeric)
					return r.included;
				if (r.pending())
					return r.draft == Disposition::Auto ? tr("Auto pending")
					       : r.desiredIncluded()        ? tr("Include pending")
					                                    : tr("Exclude pending");
				return r.included ? tr("Included") : r.replaced ? tr("Replaced") : tr("Excluded");
			default:
				return {};
		}
	}

	void reset(const orbital_core::FitPipelineResult& result)
	{
		beginResetModel();
		rows.clear();
		rows.reserve(result.observations.all.size());
		QSet<qsizetype> included;
		for (const auto i : result.fit.keptObsIndices)
			included.insert(static_cast<qsizetype>(i));
		const auto& allResiduals = result.fit.allResiduals.data();
		const bool completeResiduals = allResiduals.size() == result.observations.all.size();
		std::unordered_map<orbital_core::ObservationId, Disposition> overrides;
		for (const auto& choice : result.observationOverrides)
			overrides.emplace(choice.id, choice.disposition);
		for (std::size_t i = 0; i < result.observations.all.size(); ++i) {
			const auto& o = result.observations.all[i];
			Row r;
			r.id = o.id;
			r.radar = o.radar;
			const auto manual = overrides.find(o.id);
			if (manual != overrides.end())
				r.accepted = r.draft = manual->second;
			r.jd = o.JD;
			r.epoch = QString::fromStdString(Convertors::convertJDToGregorian(o.JD));
			r.timeScale = timeScaleName(o.timeScale);
			r.station = QString::fromStdString(o.obsCode);
			r.technique = o.radar ? tr("Radar") : techniqueName(o.note2);
			r.catalog = catalogName(o);
			r.program = !o.adesProgramCode.empty()                      ? QString::fromStdString(o.adesProgramCode)
			            : o.programCode != ' ' && o.programCode != '\0' ? QString(QChar::fromLatin1(o.programCode))
			                                                            : QString();
			r.debiasEnabled = result.observations.preparation && result.observations.preparation->debiasRequested;
			r.biasRaArcsec = o.biasRaArcsec;
			r.biasDecArcsec = o.biasDecArcsec;
			r.ra = o.submittedRaRad();
			r.dec = o.submittedDecRad();
			r.included = result.fit.keptObsIndices.empty() || included.contains(static_cast<qsizetype>(i));
			r.replaced = o.note2 == 'X' || o.note2 == 'x';
			if (completeResiduals)
				r.residual = allResiduals[i];
			rows.push_back(std::move(r));
		}
		const auto fittedCount = result.fit.keptObsIndices.empty() ? rows.size() : result.fit.keptObsIndices.size();
		if (!completeResiduals && result.fit.residuals.count() == fittedCount) {
			for (std::size_t j = 0; j < fittedCount; ++j) {
				const auto i = result.fit.keptObsIndices.empty() ? j : result.fit.keptObsIndices[j];
				if (i < rows.size())
					rows[i].residual = result.fit.residuals.data()[j];
			}
		}
		rangeRates =
		    std::ranges::any_of(rows, [](const auto& row) { return row.radar && row.radar->frequencyIsReference; });
		endResetModel();
	}

	void notifyEdits(int row)
	{
		emit dataChanged(index(row, 0), index(row, columnCount() - 1),
		                 {Qt::DisplayRole, desiredIncludedRole, pendingRole, manualRole});
	}
};

ObservationTableModel::ObservationTableModel(QObject* parent)
    : QSortFilterProxyModel(parent), m_rows(std::make_unique<ObservationRows>())
{
	setSourceModel(m_rows.get());
	setSortRole(numericRole);
	setDynamicSortFilter(true);
	connect(this, &QAbstractItemModel::modelReset, this, &ObservationTableModel::countChanged);
	connect(this, &QAbstractItemModel::rowsInserted, this, &ObservationTableModel::countChanged);
	connect(this, &QAbstractItemModel::rowsRemoved, this, &ObservationTableModel::countChanged);
	connect(this, &QAbstractItemModel::layoutChanged, this, &ObservationTableModel::selectionChanged);
	connect(this, &QAbstractItemModel::rowsInserted, this, &ObservationTableModel::selectionChanged);
	connect(this, &QAbstractItemModel::rowsRemoved, this, &ObservationTableModel::selectionChanged);
}

ObservationTableModel::~ObservationTableModel() = default;

void ObservationTableModel::setResult(const orbital_core::FitPipelineResult& result)
{
	const auto selectedId =
	    m_selectedSourceRow >= 0 ? m_rows->rows[m_selectedSourceRow].id : orbital_core::ObservationId();
	m_selectedSourceRow = -1;
	m_sourceFormat = QString::fromStdString(result.observations.format);
	m_rows->reset(result);
	if (radarMode() && radarCount() == 0)
		setRadarMode(false);
	for (int i = 0; i < m_rows->rowCount(); ++i)
		if (m_rows->rows[i].id == selectedId)
			m_selectedSourceRow = i;
	std::vector<StationResidualSample> samples;
	samples.reserve(m_rows->rows.size());
	for (const auto& row : m_rows->rows) {
		samples.push_back({row.station, row.included,
		                   row.residual && !row.radar
		                       ? std::optional(QPointF(row.residual->dRA_arcsec, row.residual->dDEC_arcsec))
		                       : std::nullopt});
	}
	m_stationResiduals.setObservations(samples);
	if (!m_stationFilter.isEmpty() && !m_stationResiduals.stations().contains(m_stationFilter))
		setStationFilter({});
	emit countChanged();
	emit selectionChanged();
	emit editsChanged();
}

int ObservationTableModel::pendingCount() const
{
	return static_cast<int>(std::ranges::count_if(m_rows->rows, &ObservationRows::Row::pending));
}

int ObservationTableModel::manualCount() const
{
	return static_cast<int>(
	    std::ranges::count_if(m_rows->rows, [](const auto& r) { return r.draft != Disposition::Auto; }));
}

void ObservationTableModel::setEditingEnabled(bool enabled)
{
	if (m_editingEnabled == enabled)
		return;
	m_editingEnabled = enabled;
	emit editingEnabledChanged();
}

std::vector<orbital_core::ObservationOverride> ObservationTableModel::overrides() const
{
	std::vector<orbital_core::ObservationOverride> result;
	for (const auto& row : m_rows->rows)
		if (row.draft != Disposition::Auto)
			result.push_back({row.id, row.draft});
	return result;
}

void ObservationTableModel::setIncluded(int sourceRow, bool included)
{
	if (!m_editingEnabled || sourceRow < 0 || sourceRow >= totalCount())
		return;
	auto& row = m_rows->rows[sourceRow];
	const auto draft = included == row.included ? row.accepted : included ? Disposition::Include : Disposition::Exclude;
	if (row.draft == draft)
		return;
	row.draft = draft;
	m_rows->notifyEdits(sourceRow);
	emit editsChanged();
	emit selectionChanged();
}

void ObservationTableModel::resetManualSelection()
{
	if (!m_editingEnabled)
		return;
	for (int i = 0; i < totalCount(); ++i) {
		m_rows->rows[i].draft = Disposition::Auto;
		m_rows->notifyEdits(i);
	}
	emit editsChanged();
	emit selectionChanged();
}

void ObservationTableModel::discardEdits()
{
	if (!m_editingEnabled)
		return;
	for (int i = 0; i < totalCount(); ++i) {
		auto& row = m_rows->rows[i];
		row.draft = row.accepted;
		m_rows->notifyEdits(i);
	}
	emit editsChanged();
	emit selectionChanged();
}

void ObservationTableModel::setFilterText(const QString& value)
{
	if (m_filter == value)
		return;
	m_filter = value;
	invalidateFilter();
	emit filterChanged();
	emit countChanged();
}

int ObservationTableModel::totalCount() const { return m_rows->rowCount(); }

int ObservationTableModel::radarCount() const
{
	return static_cast<int>(std::ranges::count_if(m_rows->rows, [](const auto& r) { return r.radar.has_value(); }));
}

bool ObservationTableModel::radarMode() const { return m_rows->radarMode; }

void ObservationTableModel::setRadarMode(bool radar)
{
	if (m_rows->radarMode == radar)
		return;
	m_rows->radarMode = radar;
	m_selectedSourceRow = -1;
	m_stationFilter.clear();
	invalidateFilter();
	emit m_rows->headerDataChanged(Qt::Horizontal, 0, m_rows->columnCount() - 1);
	emit filterChanged();
	emit countChanged();
	emit selectionChanged();
}

int ObservationTableModel::selectedProxyRow() const
{
	return mapFromSource(m_rows->index(m_selectedSourceRow, 0)).row();
}

void ObservationTableModel::clearFilters()
{
	if (m_filter.isEmpty() && m_stationFilter.isEmpty() && !m_excludedOnly)
		return;
	m_filter.clear();
	m_stationFilter.clear();
	m_excludedOnly = false;
	invalidateFilter();
	emit filterChanged();
	emit countChanged();
}

void ObservationTableModel::setExcludedOnly(bool value)
{
	if (m_excludedOnly == value)
		return;
	m_excludedOnly = value;
	invalidateFilter();
	emit filterChanged();
	emit countChanged();
}

void ObservationTableModel::setStationFilter(const QString& value)
{
	if (m_stationFilter == value)
		return;
	m_stationFilter = value;
	invalidateFilter();
	emit filterChanged();
	emit countChanged();
}

void ObservationTableModel::showStation(const QString& station)
{
	setRadarMode(false);
	m_filter.clear();
	m_excludedOnly = false;
	m_stationFilter = station;
	m_selectedSourceRow = -1;
	invalidateFilter();
	emit filterChanged();
	emit countChanged();
	emit selectionChanged();
}

bool ObservationTableModel::filterAcceptsRow(int row, const QModelIndex&) const
{
	const auto& r = m_rows->rows[row];
	return r.radar.has_value() == radarMode() && (!m_excludedOnly || !r.included) &&
	       (m_stationFilter.isEmpty() || r.station == m_stationFilter) &&
	       (m_filter.trimmed().isEmpty() || (r.epoch + " " + r.station + " " + r.program + " " + r.catalog + " " +
	                                         r.technique + " " + QString::number(row + 1))
	                                            .contains(m_filter.trimmed(), Qt::CaseInsensitive));
}

void ObservationTableModel::selectRow(int proxyRow)
{
	if (proxyRow < 0 || proxyRow >= rowCount())
		return;
	selectSourceRow(mapToSource(index(proxyRow, 0)).row());
}

void ObservationTableModel::selectSourceRow(int sourceRow)
{
	if (sourceRow < 0 || sourceRow >= m_rows->rowCount())
		return;
	setRadarMode(m_rows->rows[sourceRow].radar.has_value());
	m_selectedSourceRow = sourceRow;
	emit selectionChanged();
}

QVariantMap ObservationTableModel::rowData(int sourceRow) const
{
	if (sourceRow < 0 || sourceRow >= m_rows->rowCount())
		return {};
	const auto& r = m_rows->rows[sourceRow];
	const auto* weight = r.residual && r.residual->astrometricWeight ? &*r.residual->astrometricWeight : nullptr;
	const bool replacedAndExcluded = r.replaced && !r.included;
	const bool weighted = !replacedAndExcluded && weight && weight->weighted;
	QVariantMap data{
	    {"sourceRow", sourceRow},
	    {"number", sourceRow + 1},
	    {"epoch", r.epoch},
	    {"timeScale", r.timeScale},
	    {"radar", r.radar.has_value()},
	    {"station", r.station},
	    {"included", r.included},
	    {"replaced", r.replaced},
	    {"desiredIncluded", r.desiredIncluded()},
	    {"pending", r.pending()},
	    {"manual", r.draft != Disposition::Auto},
	    {"technique", r.technique},
	    {"catalog", r.catalog},
	    {"program", r.program},
	    {"debiasEnabled", r.debiasEnabled},
	    {"biasRa", number(r.biasRaArcsec)},
	    {"biasDec", number(r.biasDecArcsec)},
	    {"weighting", replacedAndExcluded ? tr("Replaced by MPC")
	                  : weight            ? (weighted ? tr("Weighted") : tr("Unweighted"))
	                                      : tr("Not available")},
	    {"sigmaRa", weighted ? number(weight->sigmaRaArcsec) : QStringLiteral("-")},
	    {"sigmaDec", weighted ? number(weight->sigmaDecArcsec) : QStringLiteral("-")},
	    {"chi", weighted && weight->chiSquared >= 0 ? number(std::sqrt(weight->chiSquared), 2) : QStringLiteral("-")},
	    {"postFitChi", weighted && weight->postFitChiSquared && *weight->postFitChiSquared >= 0
	                       ? number(std::sqrt(*weight->postFitChiSquared), 2)
	                       : QStringLiteral("-")},
	    {"ra", r.radar ? QStringLiteral("-") : FormatUtils::ra(r.ra)},
	    {"dec", r.radar ? QStringLiteral("-") : declination(r.dec)},
	    {"raResidual", r.residual ? number(r.residual->dRA_arcsec) : QStringLiteral("-")},
	    {"decResidual", r.residual ? number(r.residual->dDEC_arcsec) : QStringLiteral("-")}};
	if (r.radar) {
		const auto* residual = r.residual && r.residual->radar ? &*r.residual->radar : nullptr;
		data["transmitter"] = QString::fromStdString(r.radar->transmitter);
		data["frequency"] = r.radar->frequencyIsReference ? QStringLiteral("-") : number(r.radar->frequencyHz / 1e6, 6);
		data["returnType"] =
		    r.radar->returnType == orbital_core::obs::RadarReturn::CenterOfMass ? tr("Center of mass") : tr("Surface");
		QVariantList channels;
		const auto add = [&](const QString& label, const auto& measurement, const auto* channel, double scale,
		                     int precision) {
			if (!measurement)
				return;
			channels.push_back(QVariantMap{
			    {"label", label},
			    {"observed", number(measurement->value * scale, precision)},
			    {"sigma", number(measurement->sigma * std::abs(scale), precision)},
			    {"computed", channel ? number(channel->computed * scale, precision) : QStringLiteral("-")},
			    {"residual", channel ? number(channel->residual * scale, precision) : QStringLiteral("-")}});
		};
		add(tr("Delay, us"), r.radar->delaySeconds,
		    residual && residual->delaySeconds ? &*residual->delaySeconds : nullptr, 1e6, 3);
		add(r.radar->frequencyIsReference ? tr("Range rate, m/s") : tr("Doppler, Hz"), r.radar->dopplerHz,
		    residual && residual->dopplerHz ? &*residual->dopplerHz : nullptr,
		    r.radar->frequencyIsReference ? rangeRateScale(*r.radar) : 1.0, 6);
		data["radarChannels"] = channels;
		data["chi"] = residual ? number(std::sqrt(residual->chiSquared), 2) : QStringLiteral("-");
		data["postFitChi"] = residual && residual->postFitChiSquared
		                         ? number(std::sqrt(*residual->postFitChiSquared), 2)
		                         : QStringLiteral("-");
	}
	return data;
}

QVariantMap ObservationTableModel::selection() const { return rowData(m_selectedSourceRow); }

QRectF ObservationTableModel::plotBounds(bool timeMode) const
{
	double xLimit = 0.1, yLimit = 0.1;
	double start = std::numeric_limits<double>::infinity();
	double end = -start;
	for (const auto& row : m_rows->rows) {
		if (row.radar || !row.residual || row.residual->radar || !std::isfinite(row.residual->dRA_arcsec) ||
		    !std::isfinite(row.residual->dDEC_arcsec))
			continue;
		xLimit = std::max(xLimit, std::abs(row.residual->dRA_arcsec) * 1.15);
		yLimit = std::max(yLimit, std::abs(row.residual->dDEC_arcsec) * 1.15);
		if (std::isfinite(row.jd)) {
			start = std::min(start, row.jd - 2400000.5);
			end = std::max(end, row.jd - 2400000.5);
		}
	}
	if (!timeMode)
		return {-xLimit, -yLimit, 2 * xLimit, 2 * yLimit};
	if (!std::isfinite(start)) {
		start = 0;
		end = 1;
	}
	const double margin = std::max(0.01, (end - start) * 0.04);
	yLimit = std::max(xLimit, yLimit);
	return {start - margin, -yLimit, end - start + 2 * margin, 2 * yLimit};
}

void ObservationTableModel::fillSeries(QObject* object, bool included, bool timeMode, bool dec, bool selectedOnly,
                                       double timeOriginMjd) const
{
	auto* series = qobject_cast<QXYSeries*>(object);
	if (!series)
		return;
	QList<QPointF> points;
	for (int i = 0; i < m_rows->rowCount(); ++i) {
		const auto& r = m_rows->rows[i];
		if (r.radar || !r.residual || r.residual->radar ||
		    (selectedOnly ? i != m_selectedSourceRow : r.included != included))
			continue;
		const double x = timeMode ? r.jd - 2400000.5 - timeOriginMjd : r.residual->dRA_arcsec;
		const double y = timeMode && !dec ? r.residual->dRA_arcsec : r.residual->dDEC_arcsec;
		if (std::isfinite(x) && std::isfinite(y))
			points.append({x, y});
	}
	series->replace(points);
}

void ObservationTableModel::selectPoint(double x, double y, bool timeMode, bool dec)
{
	selectSourceRow(closestSourceRow(x, y, timeMode, dec));
}

QVariantMap ObservationTableModel::pointData(double x, double y, bool timeMode, bool dec, bool included) const
{
	return rowData(closestSourceRow(x, y, timeMode, dec, included));
}

int ObservationTableModel::closestSourceRow(double x, double y, bool timeMode, bool dec,
                                            std::optional<bool> included) const
{
	if (!std::isfinite(x) || !std::isfinite(y))
		return -1;
	double closest = std::numeric_limits<double>::infinity();
	int selected = -1;
	const auto bounds = plotBounds(timeMode);
	for (int i = 0; i < m_rows->rowCount(); ++i) {
		const auto& r = m_rows->rows[i];
		if (r.radar || !r.residual || r.residual->radar || (included && r.included != *included))
			continue;
		const double px = timeMode ? r.jd - 2400000.5 : r.residual->dRA_arcsec;
		const double py = timeMode && !dec ? r.residual->dRA_arcsec : r.residual->dDEC_arcsec;
		const double distance = std::hypot((px - x) / bounds.width(), (py - y) / bounds.height());
		if (distance < closest) {
			closest = distance;
			selected = i;
		}
	}
	return selected;
}

void ObservationTableModel::sortColumn(int column, bool descending)
{
	if (column < 0 || column >= columnCount())
		return;
	if (sortedColumn() == column && sortDescending() == descending)
		return;
	sort(column, descending ? Qt::DescendingOrder : Qt::AscendingOrder);
	emit sortingChanged();
}

QString ObservationTableModel::selectionTitle() const
{
	return m_selectedSourceRow < 0 ? tr("Observation") : tr("Observation %1").arg(m_selectedSourceRow + 1);
}

QString ObservationTableModel::selectionDetails() const
{
	if (m_selectedSourceRow < 0 || m_selectedSourceRow >= m_rows->rowCount())
		return {};
	const auto& r = m_rows->rows[m_selectedSourceRow];
	if (r.radar)
		return tr("%1 %2\nRadar: %3 -> %4")
		    .arg(r.epoch, r.timeScale, QString::fromStdString(r.radar->transmitter), r.station);
	QString text = tr("%1 %2\nStation: %3\n\nRA: %4\nDec: %5\n\n%6")
	                   .arg(r.epoch, r.timeScale, r.station, FormatUtils::ra(r.ra), declination(r.dec),
	                        r.included   ? tr("Included in fit")
	                        : r.replaced ? tr("Replaced by MPC")
	                                     : tr("Excluded from fit"));
	if (r.residual)
		text += tr("\n\ndRA cos Dec: %1 arcsec\ndDec: %2 arcsec")
		            .arg(number(r.residual->dRA_arcsec), number(r.residual->dDEC_arcsec));
	return text;
}
