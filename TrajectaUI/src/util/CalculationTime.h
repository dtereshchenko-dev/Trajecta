#pragma once

#include <QDateTime>
#include <QTimeZone>
#include <orbital_core/FitResult.h>
#include <orbital_core/TimeScale.h>

#include <cmath>
#include <optional>

namespace CalculationTime {

	inline QString formatUtc(double mjd)
	{
		if (!std::isfinite(mjd) || mjd < -678575 || mjd >= 2973484)
			return {};
		const auto date = QDateTime::fromMSecsSinceEpoch(
		    static_cast<qint64>(std::llround((mjd - 40587.0) * 86400000.0)), QTimeZone::UTC);
		if (!date.isValid() || date.date().year() < 1 || date.date().year() > 9999)
			return {};
		return date.toString(QStringLiteral("yyyy-MM-dd HH:mm:ss.zzz"));
	}

	inline std::optional<double> defaultStartUtc(const orbital_core::FitResult& fit)
	{
		const auto& fitOrbit = fit.fitEpochOrbit ? fit.fitEpochOrbit : fit.orbit;
		if (!fitOrbit || !std::isfinite(fitOrbit->epochMjd()))
			return std::nullopt;
		const double fitEpoch = fitOrbit->epochMjd();
		for (const auto& orbit : {fitOrbit, fit.orbit}) {
			if (!orbit || !std::isfinite(orbit->epochMjd()) || orbit->epochMjd() < fitEpoch)
				continue;
			const double utc = orbital_core::time::fromMjdTdb(orbit->epochMjd(), orbital_core::time::Scale::UTC);
			if (std::isfinite(utc) && utc >= -678575 && utc < 2973484)
				return utc;
		}
		const double now = 40587.0 + QDateTime::currentMSecsSinceEpoch() / 86400000.0;
		return orbital_core::time::toMjdTdb(now, orbital_core::time::Scale::UTC) >= fitEpoch
		           ? std::optional<double>(now)
		           : std::nullopt;
	}

}  // namespace CalculationTime
