#include "CovarianceTableModel.h"
#include <orbital_core/Constants.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
	constexpr double missing = std::numeric_limits<double>::quiet_NaN();
	constexpr int rhoRole = Qt::UserRole + 1;
}  // namespace

int CovarianceTableModel::rowCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : dimension(); }

int CovarianceTableModel::columnCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : dimension(); }

double CovarianceTableModel::rho(int row, int column) const
{
	const double a = m_matrix(row, row), b = m_matrix(column, column);
	if (!(a > 0 && b > 0))
		return missing;
	const double value = m_matrix(row, column) / std::sqrt(a) / std::sqrt(b);
	return std::isfinite(value) && std::abs(value) <= 1.0 + 1e-12 ? std::clamp(value, -1.0, 1.0) : missing;
}

QVariant CovarianceTableModel::data(const QModelIndex& index, int role) const
{
	if (!index.isValid() || index.row() < 0 || index.column() < 0 || index.row() >= dimension() ||
	    index.column() >= dimension())
		return {};
	if (role != Qt::DisplayRole && role != Qt::UserRole && role != rhoRole)
		return {};
	const double value =
	    m_correlation || role == rhoRole ? rho(index.row(), index.column()) : m_matrix(index.row(), index.column());
	if (!std::isfinite(value))
		return {};
	return role == Qt::DisplayRole ? QVariant(QString::number(value, m_correlation ? 'f' : 'e', m_correlation ? 3 : 1))
	                               : QVariant(value);
}

QVariant CovarianceTableModel::headerData(int section, Qt::Orientation, int role) const
{
	return role == Qt::DisplayRole ? QVariant(m_labels.value(section)) : QVariant();
}

QHash<int, QByteArray> CovarianceTableModel::roleNames() const
{
	return {{Qt::DisplayRole, "display"}, {Qt::UserRole, "numeric"}, {rhoRole, "rho"}};
}

const CovarianceTableModel::Diagnostics* CovarianceTableModel::diagnostics() const
{
	const auto& view = m_equinoctial ? m_equ : m_car;
	return view ? &*view : nullptr;
}

void CovarianceTableModel::setResult(const orbital_core::FitResult& result)
{
	using namespace orbital_core::uncertainty;
	m_car = describeCovariance(result, CovarianceBasis::Cartesian);
	m_equ = describeCovariance(result, CovarianceBasis::Equinoctial);
	m_equinoctial = m_equ.has_value();
	m_nonGravLabels.clear();
	m_nonGravUnits.clear();
	const QStringList labels{"A1", "A2", "A3", "DT"};
	for (int i = 0; i < result.finalNonGrav.numFitParams(); ++i) {
		const int parameter = static_cast<int>(result.finalNonGrav.fitParameterAt(i));
		m_nonGravLabels.append(labels.value(parameter));
		m_nonGravUnits.append(parameter == 3 ? QStringLiteral("days") : QStringLiteral("AU/day\u00b2"));
	}
	m_selectedRow = m_selectedColumn = -1;
	rebuild();
}

void CovarianceTableModel::rebuild()
{
	beginResetModel();
	m_matrix.resize(0, 0);
	m_labels.clear();
	m_units.clear();
	if (const auto* view = diagnostics()) {
		m_matrix = view->covariance;
		if (m_equinoctial) {
			m_labels = {"a", "h", "k", "p", "q", QStringLiteral("\u03bb")};
			m_units = {"AU", "", "", "", "", "deg"};
			m_matrix.row(5) *= Constants::RAD_TO_DEG;
			m_matrix.col(5) *= Constants::RAD_TO_DEG;
		} else {
			m_labels = {"X", "Y", "Z", "Vx", "Vy", "Vz"};
			m_units = {"AU", "AU", "AU", "AU/day", "AU/day", "AU/day"};
		}
		m_labels.append(m_nonGravLabels);
		m_units.append(m_nonGravUnits);
	}
	m_selectedRow = dimension() ? std::clamp(m_selectedRow, 0, dimension() - 1) : -1;
	m_selectedColumn = dimension() ? std::clamp(m_selectedColumn, 0, dimension() - 1) : -1;
	endResetModel();
	emit changed();
	emit selectionChanged();
}

void CovarianceTableModel::setCorrelation(bool value)
{
	if (value == m_correlation)
		return;
	m_correlation = value;
	if (dimension())
		emit dataChanged(index(0, 0), index(dimension() - 1, dimension() - 1));
	emit changed();
}

void CovarianceTableModel::setEquinoctial(bool value)
{
	if (value == m_equinoctial || (value && !m_equ))
		return;
	m_equinoctial = value;
	rebuild();
}

QString CovarianceTableModel::epoch() const
{
	const auto* view = diagnostics();
	return view ? QString::number(view->epochMjdTdb, 'f', 9) + QStringLiteral(" MJD TDB") : QString();
}

double CovarianceTableModel::objectiveScale() const
{
	const auto* view = diagnostics();
	return view ? view->objectiveScale : missing;
}

void CovarianceTableModel::selectCell(int row, int column)
{
	if (row < 0 || column < 0 || row >= dimension() || column >= dimension() ||
	    (row == m_selectedRow && column == m_selectedColumn))
		return;
	m_selectedRow = row;
	m_selectedColumn = column;
	emit selectionChanged();
}

QVariantMap CovarianceTableModel::selection() const
{
	if (m_selectedRow < 0 || m_selectedColumn < 0)
		return {};
	const int r = m_selectedRow, c = m_selectedColumn;
	const auto sigma = [&](int i) { return m_matrix(i, i) >= 0 ? std::sqrt(m_matrix(i, i)) : missing; };
	return {{"rowLabel", m_labels[r]},  {"columnLabel", m_labels[c]}, {"rowUnit", m_units[r]},
	        {"columnUnit", m_units[c]}, {"correlation", rho(r, c)},   {"covariance", m_matrix(r, c)},
	        {"rowSigma", sigma(r)},     {"columnSigma", sigma(c)}};
}

double CovarianceTableModel::weakSigma() const
{
	const auto* view = diagnostics();
	if (!view || !view->spectrum || view->spectrum->signedSigmas.size() == 0)
		return missing;
	const auto& sigmas = view->spectrum->signedSigmas;
	const double last = sigmas(sigmas.size() - 1);
	if (!(last > 0) || sigmas(0) < -last * 1e-6)
		return missing;
	return last;
}

QVariantList CovarianceTableModel::weakDirection() const
{
	if (!std::isfinite(weakSigma()))
		return {};
	QVariantList rows;
	const auto& axis = diagnostics()->spectrum->weakDirection;
	for (int i = 0; i < dimension(); ++i)
		rows.append(QVariant::fromValue(QVariantList{m_labels[i], axis(i)}));
	return rows;
}

QString CovarianceTableModel::basisUnits() const
{
	QString units = m_equinoctial ? QStringLiteral("EQU: a [AU], h/k/p/q [1], \u03bb [deg]")
	                              : QStringLiteral("CAR: X/Y/Z [AU], Vx/Vy/Vz [AU/day]");
	for (int i = 0; i < m_nonGravLabels.size(); ++i)
		units += QStringLiteral("; %1 [%2]").arg(m_nonGravLabels[i], m_nonGravUnits[i]);
	return units;
}
