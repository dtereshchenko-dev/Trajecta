#pragma once
#include <QAbstractTableModel>
#include <QQmlEngine>
#include <orbital_core/hazard/IntervalMonteCarlo.h>

class RiskEventTableModel : public QAbstractTableModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the risk monitor")
	Q_PROPERTY(int count READ count NOTIFY countChanged)
   public:
	explicit RiskEventTableModel(QObject* parent = nullptr) : QAbstractTableModel(parent) {}

	int count() const { return static_cast<int>(m_rows.size()); }

	int rowCount(const QModelIndex& parent = {}) const override { return parent.isValid() ? 0 : count(); }

	int columnCount(const QModelIndex& parent = {}) const override { return parent.isValid() ? 0 : 7; }

	QVariant data(const QModelIndex&, int role = Qt::DisplayRole) const override;
	QVariant headerData(int, Qt::Orientation, int role = Qt::DisplayRole) const override;
	QHash<int, QByteArray> roleNames() const override;
	void setEvents(const std::vector<orbital_core::hazard::McEvent>& events);

   private:
	struct Row {
		QStringList cells;
		QString detail;
		bool entry{};
	};

	std::vector<Row> m_rows;
   signals:
	void countChanged();
};
