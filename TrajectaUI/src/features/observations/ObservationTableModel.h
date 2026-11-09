#pragma once

#include <QSortFilterProxyModel>
#include <QQmlEngine>
#include <QRectF>
#include <QVariantMap>
#include <memory>
#include <optional>
#include <orbital_core/FitPipeline.h>
#include "StationResidualModel.h"

class ObservationRows;

class ObservationTableModel : public QSortFilterProxyModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the solution")
	Q_PROPERTY(QString filterText READ filterText WRITE setFilterText NOTIFY filterChanged)
	Q_PROPERTY(QString stationFilter READ stationFilter WRITE setStationFilter NOTIFY filterChanged)
	Q_PROPERTY(StationResidualModel* stationResiduals READ stationResiduals CONSTANT)
	Q_PROPERTY(bool excludedOnly READ excludedOnly WRITE setExcludedOnly NOTIFY filterChanged)
	Q_PROPERTY(int count READ count NOTIFY countChanged)
	Q_PROPERTY(int totalCount READ totalCount NOTIFY countChanged)
	Q_PROPERTY(int radarCount READ radarCount NOTIFY countChanged)
	Q_PROPERTY(bool radarMode READ radarMode WRITE setRadarMode NOTIFY filterChanged)
	Q_PROPERTY(QString sourceFormat READ sourceFormat NOTIFY countChanged)
	Q_PROPERTY(int selectedSourceRow READ selectedSourceRow NOTIFY selectionChanged)
	Q_PROPERTY(int selectedProxyRow READ selectedProxyRow NOTIFY selectionChanged)
	Q_PROPERTY(int sortedColumn READ sortedColumn NOTIFY sortingChanged)
	Q_PROPERTY(bool sortDescending READ sortDescending NOTIFY sortingChanged)
	Q_PROPERTY(QString selectionTitle READ selectionTitle NOTIFY selectionChanged)
	Q_PROPERTY(QString selectionDetails READ selectionDetails NOTIFY selectionChanged)
	Q_PROPERTY(QVariantMap selection READ selection NOTIFY selectionChanged)
	Q_PROPERTY(int pendingCount READ pendingCount NOTIFY editsChanged)
	Q_PROPERTY(int manualCount READ manualCount NOTIFY editsChanged)
	Q_PROPERTY(bool editingEnabled READ editingEnabled WRITE setEditingEnabled NOTIFY editingEnabledChanged)

   public:
	explicit ObservationTableModel(QObject* parent = nullptr);
	~ObservationTableModel() override;
	void setResult(const orbital_core::FitPipelineResult& result);

	QString filterText() const { return m_filter; }

	QString stationFilter() const { return m_stationFilter; }

	StationResidualModel* stationResiduals() { return &m_stationResiduals; }

	bool excludedOnly() const { return m_excludedOnly; }

	int count() const { return rowCount(); }

	int totalCount() const;
	int radarCount() const;
	bool radarMode() const;
	void setRadarMode(bool radar);

	QString sourceFormat() const { return m_sourceFormat; }

	int selectedSourceRow() const { return m_selectedSourceRow; }

	int selectedProxyRow() const;

	int sortedColumn() const { return QSortFilterProxyModel::sortColumn(); }

	bool sortDescending() const { return sortOrder() == Qt::DescendingOrder; }

	QString selectionTitle() const;
	QString selectionDetails() const;
	QVariantMap selection() const;
	int pendingCount() const;
	int manualCount() const;

	bool editingEnabled() const { return m_editingEnabled; }

	void setEditingEnabled(bool enabled);
	std::vector<orbital_core::ObservationOverride> overrides() const;
	Q_INVOKABLE void setIncluded(int sourceRow, bool included);
	Q_INVOKABLE void resetManualSelection();
	Q_INVOKABLE void discardEdits();
	void setFilterText(const QString& value);
	void setStationFilter(const QString& value);
	Q_INVOKABLE void showStation(const QString& station);
	void setExcludedOnly(bool value);
	Q_INVOKABLE void selectRow(int proxyRow);
	Q_INVOKABLE void sortColumn(int column, bool descending);
	Q_INVOKABLE void selectSourceRow(int sourceRow);
	Q_INVOKABLE void clearFilters();
	Q_INVOKABLE QRectF plotBounds(bool timeMode) const;
	Q_INVOKABLE void fillSeries(QObject* series, bool included, bool timeMode, bool declination,
	                            bool selectedOnly = false, double timeOriginMjd = 0) const;
	Q_INVOKABLE void selectPoint(double x, double y, bool timeMode, bool declination);
	Q_INVOKABLE QVariantMap pointData(double x, double y, bool timeMode, bool declination, bool included) const;

   signals:
	void filterChanged();
	void countChanged();
	void selectionChanged();
	void sortingChanged();
	void editsChanged();
	void editingEnabledChanged();

   protected:
	bool filterAcceptsRow(int row, const QModelIndex& parent) const override;

   private:
	void beginObservationFilterChange();
	void endObservationFilterChange();
	QVariantMap rowData(int sourceRow) const;
	int closestSourceRow(double x, double y, bool timeMode, bool declination, std::optional<bool> included = {}) const;
	std::unique_ptr<ObservationRows> m_rows;
	StationResidualModel m_stationResiduals;
	QString m_filter;
	QString m_stationFilter;
	QString m_sourceFormat;
	bool m_excludedOnly{};
	bool m_editingEnabled{true};
	int m_selectedSourceRow{-1};
};
