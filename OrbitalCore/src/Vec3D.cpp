#include "orbital_core/Vec3D.h"
#include <string>

#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>

std::ostream& operator<<(std::ostream& os, const Vec3D& vec)
{
	os << "(" << std::fixed << std::setprecision(10) << vec.x << ", " << vec.y << ", " << vec.z << ")";
	os.unsetf(std::ios::fixed);
	return os;
}

std::string Vec3D::print() const
{
	std::ostringstream oss;
	oss << "(" << std::fixed << std::setprecision(10) << x << ", " << y << ", " << z << ")";
	oss.unsetf(std::ios::fixed);
	return oss.str();
}
