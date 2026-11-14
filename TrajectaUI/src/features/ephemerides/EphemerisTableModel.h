#pragma once

#include <QAbstractTableModel>
#include <QQmlEngine>
#include <QVariantList>
#include <orbital_core/Ephemeris.h>

#include <array>
#include <utility>
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
	Q_PROPERTY(QVariantList columnOptions READ columnOptions NOTIFY changed)

   public:
	enum class Layout { Spherical, Cartesian };

	enum Column {
		Date,
		Longitude,
		Latitude,
		Distance,
		SolarDistance,
		Magnitude,
		Elongation,
		LongitudeRate,
		LatitudeRate,
		DistanceRate,
		SolarDistanceRate,
		LightTime,
		Phase,
		SkyMotion,
		SkyPositionAngle,
		LongitudeUncertainty,
		LatitudeUncertainty,
		EllipseMajor,
		EllipseMinor,
		EllipseAngle,
		EllipseArea,
		SkyUncertainty,
		RangeUncertainty,
		RangeRateUncertainty,
		X,
		Y,
		Z,
		Vx,
		Vy,
		Vz,
		ColumnCount
	};

	explicit EphemerisTableModel(QObject* parent = nullptr);

	void setRows(std::vector<orbital_core::EphemerisRow> rows, orbital_core::EphemerisPlane plane,
	             orbital_core::time::Scale scale, bool heliocentric,
	             orbital_core::EphemerisMagnitudeType magnitudeType = orbital_core::EphemerisMagnitudeType::Asteroid);
	void setLayout(Layout layout);
	void clear();

	[[nodiscard]] int count() const { return static_cast<int>(m_rows.size()); }

	[[nodiscard]] Layout layout() const { return m_layout; }

	[[nodiscard]] const std::vector<orbital_core::EphemerisRow>& rows() const { return m_rows; }

	std::vector<orbital_core::EphemerisRow> takeRows();

	[[nodiscard]] QString cell(int row, int column) const;
	[[nodiscard]] QString header(int column) const;
	[[nodiscard]] QString cellForColumn(int row, Column column) const;
	[[nodiscard]] QString headerForColumn(Column column) const;

	[[nodiscard]] int columns() const;
	[[nodiscard]] int logicalColumn(int column) const;
	[[nodiscard]] bool includesColumn(Column column) const;
	[[nodiscard]] QVariantList columnOptions() const;
	[[nodiscard]] bool uncertaintySelected() const;
	Q_INVOKABLE void setColumnEnabled(int column, bool enabled);
	Q_INVOKABLE int columnWidth(int column) const;

	int rowCount(const QModelIndex& parent = {}) const override;
	int columnCount(const QModelIndex& parent = {}) const override;
	QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

	// Calendar date of an MJD in its own scale: "yyyy-MM-dd HH:mm", with seconds
	// (and milliseconds when needed) only if some epoch requires them.
	[[nodiscard]] static QString calendar(double mjd, int secondDecimals);
	[[nodiscard]] static double skyMotion(const orbital_core::EphemerisRow& row);
	[[nodiscard]] static double skyPositionAngle(const orbital_core::EphemerisRow& row);

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
	std::array<std::array<bool, ColumnCount>, 2> m_enabled{};
};
