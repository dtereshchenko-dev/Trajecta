#pragma once

#include <QAbstractTableModel>
#include <QQmlEngine>
#include <orbital_core/hazard/CloseApproachSearch.h>

class CloseApproachTableModel : public QAbstractTableModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by CloseApproachesViewModel")
	Q_PROPERTY(int count READ count NOTIFY changed)
	Q_PROPERTY(int distanceUnit READ distanceUnit WRITE setDistanceUnit NOTIFY changed)
	Q_PROPERTY(int velocityUnit READ velocityUnit WRITE setVelocityUnit NOTIFY changed)
	Q_PROPERTY(int sortColumn READ sortColumn NOTIFY changed)
	Q_PROPERTY(bool descending READ descending NOTIFY changed)
   public:
	explicit CloseApproachTableModel(QObject* parent = nullptr) : QAbstractTableModel(parent) {}

	int count() const { return static_cast<int>(m_rows.size()); }

	int rowCount(const QModelIndex& parent = {}) const override { return parent.isValid() ? 0 : count(); }

	int columnCount(const QModelIndex& parent = {}) const override { return parent.isValid() ? 0 : 7; }

	QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
	QHash<int, QByteArray> roleNames() const override;
	void setRows(std::vector<orbital_core::hazard::NominalCloseApproach> rows);

	int distanceUnit() const { return m_distanceUnit; }

	int velocityUnit() const { return m_velocityUnit; }

	void setDistanceUnit(int value);
	void setVelocityUnit(int value);

	int sortColumn() const { return m_sortColumn; }

	bool descending() const { return m_descending; }

	Q_INVOKABLE void toggleSort(int column);
	QString cell(int row, int column) const;
	QString header(int column) const;
	QString toPlainText() const;
   signals:
	void changed();

   private:
	void sortRows();
	void reformat();
	std::vector<orbital_core::hazard::NominalCloseApproach> m_rows;
	int m_distanceUnit{};
	int m_velocityUnit{};
	int m_sortColumn{};
	bool m_descending{};
};
