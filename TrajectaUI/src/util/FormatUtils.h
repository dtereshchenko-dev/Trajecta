#pragma once

#include <QString>
#include <QVector3D>
#include <cmath>
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"

namespace FormatUtils {

	inline QString ra(double rad)
	{
		if (!std::isfinite(rad))
			return QStringLiteral("-");
		const double wrapped = std::fmod(std::fmod(rad, Constants::TWOPI) + Constants::TWOPI, Constants::TWOPI);
		const auto ticks = std::llround(wrapped * (12.0 / Constants::PI) * 360000.0) % 8640000;
		const int h = static_cast<int>(ticks / 360000);
		const int m = static_cast<int>(ticks / 6000 % 60);
		const double s = static_cast<double>(ticks % 6000) / 100.0;
		return QString("%1:%2:%3").arg(h, 2, 10, QChar(' ')).arg(m, 2, 10, QChar('0')).arg(s, 5, 'f', 2, QChar('0'));
	}

	inline QString dec(double rad)
	{
		if (!std::isfinite(rad) || std::abs(rad) > Constants::PI / 2.0)
			return QStringLiteral("-");
		const auto ticks = std::llround(std::abs(rad) * Constants::RAD_TO_DEG * 36000.0);
		const int d = static_cast<int>(ticks / 36000);
		const int m = static_cast<int>(ticks / 600 % 60);
		const double s = static_cast<double>(ticks % 600) / 10.0;
		return QString("%1%2 %3 %4")
		    .arg(std::signbit(rad) ? "-" : "+")
		    .arg(d)
		    .arg(m, 2, 10, QChar('0'))
		    .arg(s, 4, 'f', 1, QChar('0'));
	}

	/// Ecliptic (x,y,z) → Qt3D convention (x, z, -y)
	inline QVector3D eclToQt3D(const Vec3D& v)
	{
		return {static_cast<float>(v.x), static_cast<float>(v.z), static_cast<float>(-v.y)};
	}

}  // namespace FormatUtils
