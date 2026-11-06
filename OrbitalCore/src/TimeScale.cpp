#include "orbital_core/TimeScale.h"

#include "orbital_core/Types.h"
#include "earth_orientation/EarthOrientation.h"

namespace orbital_core::time {
	namespace {

		double ttMinusUtcDays(double mjdUtc) noexcept
		{
			return Convertors::utcToTdtSeconds(jdFromMjd(mjdUtc)) / 86400.0;
		}

	}  // namespace

	double toMjdTt(double mjd, Scale scale) noexcept
	{
		switch (scale) {
			case Scale::UT1:
				return mjd + ttOffsetSeconds(mjd, scale) / 86400.0;
			case Scale::UTC:
				return mjd + ttMinusUtcDays(mjd);
			case Scale::TDB: {
				double tt = mjd - Calculations::tdbMinusTtSeconds(mjd) / 86400.0;
				tt = mjd - Calculations::tdbMinusTtSeconds(tt) / 86400.0;
				return tt;
			}
			case Scale::TT:
			default:
				return mjd;
		}
	}

	double fromMjdTt(double mjdTt, Scale scale) noexcept
	{
		switch (scale) {
			case Scale::UT1: {
				double ut1;
				return earthOrientation().ttToUt1Mjd(mjdTt, ut1) ? ut1 : std::numeric_limits<double>::quiet_NaN();
			}
			case Scale::UTC: {
				double utc = mjdTt - ttMinusUtcDays(mjdTt);
				for (int i = 0; i < 2; ++i)
					utc = mjdTt - ttMinusUtcDays(utc);
				return utc;
			}
			case Scale::TDB:
				return mjdTt + Calculations::tdbMinusTtSeconds(mjdTt) / 86400.0;
			case Scale::TT:
			default:
				return mjdTt;
		}
	}

	double ttOffsetSeconds(double mjd, Scale scale) noexcept
	{
		switch (scale) {
			case Scale::UTC:
				return Convertors::utcToTdtSeconds(jdFromMjd(mjd));
			case Scale::UT1: {
				double seconds;
				return earthOrientation().ut1ToTtSeconds(mjd, seconds) ? seconds
				                                                       : std::numeric_limits<double>::quiet_NaN();
			}
			case Scale::TT:
				return 0.0;
			case Scale::TDB:
				return -Calculations::tdbMinusTtSeconds(toMjdTt(mjd, scale));
		}
		return std::numeric_limits<double>::quiet_NaN();
	}

	double toMjdTdb(double mjd, Scale scale) noexcept
	{
		return scale == Scale::TDB ? mjd : fromMjdTt(toMjdTt(mjd, scale), Scale::TDB);
	}

	double fromMjdTdb(double mjdTdb, Scale scale) noexcept
	{
		return scale == Scale::TDB ? mjdTdb : fromMjdTt(toMjdTt(mjdTdb, Scale::TDB), scale);
	}

}  // namespace orbital_core::time
