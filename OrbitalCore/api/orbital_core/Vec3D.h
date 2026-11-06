#ifndef VEC3D_H
#define VEC3D_H

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>
#include <string>

class Vec3D {
   public:
	double x;
	double y;
	double z;

	Vec3D() : x(0.0), y(0.0), z(0.0) {}

	Vec3D(const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}

	bool isNull() const { return x == 0.0 && y == 0.0 && z == 0.0; }

	double magnitude() const { return std::sqrt(x * x + y * y + z * z); }

	void nullify()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	Vec3D normalized() const
	{
		double mag = magnitude();
		return (mag > 0) ? Vec3D(x / mag, y / mag, z / mag) : Vec3D(0, 0, 0);
	}

	Vec3D operator+(const Vec3D& other) const { return Vec3D(x + other.x, y + other.y, z + other.z); }

	Vec3D operator-(const Vec3D& other) const { return Vec3D(x - other.x, y - other.y, z - other.z); }

	Vec3D operator*(double scalar) const { return Vec3D(x * scalar, y * scalar, z * scalar); }

	Vec3D operator*(const Vec3D& other) const { return Vec3D(x * other.x, y * other.y, z * other.z); }

	Vec3D operator/(double scalar) const { return Vec3D(x / scalar, y / scalar, z / scalar); }

	Vec3D& operator+=(const Vec3D& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	Vec3D& operator-=(const Vec3D& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	Vec3D& operator*=(double scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	Vec3D& operator/=(double scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		return *this;
	}

	double dot(const Vec3D& other) const { return x * other.x + y * other.y + z * other.z; }

	Vec3D cross(const Vec3D& other) const
	{
		return Vec3D(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
	}

	double& operator[](size_t index)
	{
		if (index == 0)
			return x;
		if (index == 1)
			return y;
		if (index == 2)
			return z;
		throw std::out_of_range("Index out of range for Vec3D");
	}

	const double& operator[](size_t index) const
	{
		if (index == 0)
			return x;
		if (index == 1)
			return y;
		if (index == 2)
			return z;
		throw std::out_of_range("Index out of range for Vec3D");
	}

	friend std::ostream& operator<<(std::ostream& os, const Vec3D& vec);
	std::string print() const;
};

#endif  // VEC3D_H
