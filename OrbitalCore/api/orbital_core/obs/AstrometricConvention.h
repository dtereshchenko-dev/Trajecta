#pragma once

namespace orbital_core::obs {
	// CatalogReferenced subtracts the deflection of the reference stars;
	// Absolute retains the target deflection. This is independent of star catalogue.
	enum class AstrometricConvention { CatalogReferenced, Absolute, GaiaFpr };

}  // namespace orbital_core::obs
