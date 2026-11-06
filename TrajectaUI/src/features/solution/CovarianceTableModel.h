#pragma once

#include <QAbstractTableModel>
#include <QQmlEngine>
#include <QVariantList>
#include <QVariantMap>
#include <orbital_core/uncertainty/CovarianceDiagnostics.h>

class CovarianceTableModel : public QAbstractTableModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the solution")
	Q_PROPERTY(bool correlation READ correlation WRITE setCorrelation NOTIFY changed)
	Q_PROPERTY(bool equinoctial READ equinoctial WRITE setEquinoctial NOTIFY changed)
	Q_PROPERTY(bool equinoctialAvailable READ equinoctialAvailable NOTIFY changed)
	Q_PROPERTY(int dimension READ dimension NOTIFY changed)
	Q_PROPERTY(QStringList labels READ labels NOTIFY changed)
	Q_PROPERTY(QString epoch READ epoch NOTIFY changed)
	Q_PROPERTY(double objectiveScale READ objectiveScale NOTIFY changed)
	Q_PROPERTY(int selectedRow READ selectedRow NOTIFY selectionChanged)
	Q_PROPERTY(int selectedColumn READ selectedColumn NOTIFY selectionChanged)
	Q_PROPERTY(QVariantMap selection READ selection NOTIFY selectionChanged)
	Q_PROPERTY(QVariantList weakDirection READ weakDirection NOTIFY changed)
	Q_PROPERTY(double weakSigma READ weakSigma NOTIFY changed)
	Q_PROPERTY(QString basisUnits READ basisUnits NOTIFY changed)

   public:
	using QAbstractTableModel::QAbstractTableModel;
	int rowCount(const QModelIndex& parent = {}) const override;
	int columnCount(const QModelIndex& parent = {}) const override;
	QVariant data(const QModelIndex& index, int role) const override;
	QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
	QHash<int, QByteArray> roleNames() const override;
	void setResult(const orbital_core::FitResult& result);

	bool correlation() const { return m_correlation; }

	void setCorrelation(bool value);

	bool equinoctial() const { return m_equinoctial; }

	void setEquinoctial(bool value);

	bool equinoctialAvailable() const { return m_equ.has_value(); }

	int dimension() const { return static_cast<int>(m_matrix.rows()); }

	QStringList labels() const { return m_labels; }

	QString epoch() const;
	double objectiveScale() const;

	int selectedRow() const { return m_selectedRow; }

	int selectedColumn() const { return m_selectedColumn; }

	QVariantMap selection() const;
	QVariantList weakDirection() const;
	double weakSigma() const;
	QString basisUnits() const;
	Q_INVOKABLE void selectCell(int row, int column);

   signals:
	void changed();
	void selectionChanged();

   private:
	using Diagnostics = orbital_core::uncertainty::CovarianceDiagnostics;
	const Diagnostics* diagnostics() const;
	void rebuild();
	double rho(int row, int column) const;
	orbital_core::FitParameterCovariance m_matrix;
	std::optional<Diagnostics> m_car;
	std::optional<Diagnostics> m_equ;
	QStringList m_labels;
	QStringList m_units;
	QStringList m_nonGravLabels;
	QStringList m_nonGravUnits;
	bool m_correlation{true};
	bool m_equinoctial{true};
	int m_selectedRow{-1};
	int m_selectedColumn{-1};
};
