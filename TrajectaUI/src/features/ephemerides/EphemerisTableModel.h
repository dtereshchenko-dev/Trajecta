#pragma once

#include <QAbstractTableModel>
#include <QQmlEngine>
#include <orbital_core/Ephemeris.h>

#include <vector>

// Rows of one computed ephemeris, formatted as on the Ephemerides page.
// Spherical: date, RA h m s (or ecliptic longitude), Dec (or latitude) d m s,
// Delta and r in AU, predicted V, solar elongation in degrees.
// Cartesian: date, X, Y, Z in AU and Vx, Vy, Vz in AU/day relative to the centre.
// Switching between the two layouts reformats the same rows without recomputing.
class EphemerisTableModel : public QAbstractTableModel {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by EphemerisViewModel")
	Q_PROPERTY(int count READ count NOTIFY changed)

   public:
	enum class Layout { Spherical, Cartesian };

	explicit EphemerisTableModel(QObject* parent = nullptr);

	void setRows(std::vector<orbital_core::EphemerisRow> rows, orbital_core::EphemerisPlane plane,
	             orbital_core::time::Scale scale, bool heliocentric,
	             orbital_core::EphemerisMagnitudeType magnitudeType = orbital_core::EphemerisMagnitudeType::Asteroid);
	void setLayout(Layout layout);
	void clear();

	[[nodiscard]] int count() const { return static_cast<int>(m_rows.size()); }

	[[nodiscard]] Layout layout() const { return m_layout; }

	[[nodiscard]] const std::vector<orbital_core::EphemerisRow>& rows() const { return m_rows; }

	[[nodiscard]] QString cell(int row, int column) const;
	[[nodiscard]] QString header(int column) const;

	[[nodiscard]] int columns() const { return 7; }

	int rowCount(const QModelIndex& parent = {}) const override;
	int columnCount(const QModelIndex& parent = {}) const override;
	QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

	// Calendar date of an MJD in its own scale: "yyyy-MM-dd HH:mm", with seconds
	// (and milliseconds when needed) only if some epoch requires them.
	[[nodiscard]] static QString calendar(double mjd, int secondDecimals);

   signals:
	void changed();

   private:
	std::vector<orbital_core::EphemerisRow> m_rows;
	Layout m_layout{Layout::Spherical};
	orbital_core::EphemerisPlane m_plane{orbital_core::EphemerisPlane::Equatorial};
	orbital_core::time::Scale m_scale{orbital_core::time::Scale::UTC};
	bool m_heliocentric{false};
	orbital_core::EphemerisMagnitudeType m_magnitudeType{orbital_core::EphemerisMagnitudeType::Asteroid};
	int m_secondDecimals{-1};
};
