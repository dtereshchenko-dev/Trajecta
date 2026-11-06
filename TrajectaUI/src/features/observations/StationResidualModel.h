#pragma once

#include <QAbstractListModel>
#include <QPointF>
#include <QQmlEngine>
#include <optional>
#include <span>
#include <vector>

struct StationResidualSample {
	QString station;
	bool included{};
	std::optional<QPointF> residual;
};

class StationResidualModel : public QAbstractListModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the observations")
	Q_PROPERTY(int count READ count NOTIFY changed)
	Q_PROPERTY(int included READ included NOTIFY changed)
	Q_PROPERTY(int excluded READ excluded NOTIFY changed)
	Q_PROPERTY(int missing READ missing NOTIFY changed)
	Q_PROPERTY(double meanLimit READ meanLimit NOTIFY changed)
	Q_PROPERTY(QStringList stations READ stations NOTIFY changed)
	Q_PROPERTY(SortMode sortMode READ sortMode WRITE setSortMode NOTIFY sortModeChanged)

   public:
	enum Role {
		StationRole = Qt::UserRole + 1,
		IncludedRole,
		ExcludedRole,
		ValidRole,
		RmsRole,
		MeanRaRole,
		MeanDecRole
	};

	enum SortMode { IncludedCount, Rms, MeanOffset, StationCode };
	Q_ENUM(SortMode)

	using QAbstractListModel::QAbstractListModel;
	int rowCount(const QModelIndex& parent = {}) const override;
	QVariant data(const QModelIndex& index, int role) const override;
	QHash<int, QByteArray> roleNames() const override;
	void setObservations(std::span<const StationResidualSample> observations);

	int count() const { return static_cast<int>(m_rows.size()); }

	int included() const { return m_included; }

	int excluded() const { return m_excluded; }

	int missing() const { return m_missing; }

	double meanLimit() const { return m_meanLimit; }

	QStringList stations() const { return m_stations; }

	SortMode sortMode() const { return m_sortMode; }

	void setSortMode(SortMode mode);

   signals:
	void changed();
	void sortModeChanged();

   private:
	struct Row {
		QString station;
		int included{};
		int excluded{};
		int valid{};
		double rms{};
		double meanRa{};
		double meanDec{};
	};

	void sortRows();
	std::vector<Row> m_rows;
	QStringList m_stations;
	int m_included{};
	int m_excluded{};
	int m_missing{};
	double m_meanLimit{0.1};
	SortMode m_sortMode{IncludedCount};
};
