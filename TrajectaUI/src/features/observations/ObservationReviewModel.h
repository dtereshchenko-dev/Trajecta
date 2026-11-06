#pragma once

#include <QAbstractListModel>
#include <QQmlEngine>
#include <orbital_core/ObservationSet.h>
#include <orbital_core/FitPipeline.h>

class ObservationReviewModel : public QAbstractListModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the session")
	Q_PROPERTY(int filter READ filter WRITE setFilter NOTIFY filterChanged)
	Q_PROPERTY(int positionFilter READ positionFilter CONSTANT)
	Q_PROPERTY(int differentPositionFilter READ differentPositionFilter CONSTANT)
	Q_PROPERTY(int belowHorizonFilter READ belowHorizonFilter CONSTANT)
	Q_PROPERTY(int daylightFilter READ daylightFilter CONSTANT)
	Q_PROPERTY(int totalCount READ totalCount NOTIFY changed)
	Q_PROPERTY(int shownCount READ shownCount NOTIFY filterChanged)
	Q_PROPERTY(int excludedCount READ excludedCount NOTIFY changed)
	Q_PROPERTY(int remainingOpticalCount READ remainingOpticalCount NOTIFY changed)

   public:
	explicit ObservationReviewModel(QObject* parent = nullptr);
	void setObservations(const orbital_core::ObservationSet& observations);
	void clear();
	int rowCount(const QModelIndex& parent = {}) const override;
	QVariant data(const QModelIndex& index, int role) const override;
	QHash<int, QByteArray> roleNames() const override;

	int filter() const { return m_filter; }

	void setFilter(int value);

	int positionFilter() const { return static_cast<int>(orbital_core::ObservationIssue::CoincidentPosition) + 1; }

	int differentPositionFilter() const
	{
		return static_cast<int>(orbital_core::ObservationIssue::NearbyObservationConflict) + 1;
	}

	int belowHorizonFilter() const { return static_cast<int>(orbital_core::ObservationIssue::BelowHorizon) + 1; }

	int daylightFilter() const { return static_cast<int>(orbital_core::ObservationIssue::Daylight) + 1; }

	int totalCount() const { return static_cast<int>(m_rows.size()); }

	int shownCount() const { return static_cast<int>(m_visible.size()); }

	int excludedCount() const;

	int remainingOpticalCount() const { return m_initialOpticalCount - excludedCount(); }

	Q_INVOKABLE void setIncluded(int visibleRow, bool included);
	Q_INVOKABLE void setAllShownIncluded(bool included);
	std::vector<orbital_core::ObservationOverride> overrides() const;

   signals:
	void filterChanged();
	void changed();

   private:
	struct Row {
		orbital_core::ObservationId id;
		QString station, time, coordinates, reason;
		std::vector<orbital_core::ObservationIssue> issues;
		int group{};
		bool included{true};
		bool replaced{};
	};

	void rebuildFilter();
	std::vector<Row> m_rows;
	std::vector<std::size_t> m_visible;
	int m_filter{};
	int m_initialOpticalCount{};
};
