#include "orbital_core/ObsNode.h"
#include <limits>
#include <iomanip>
#include <string>
#include <cmath>
#include "orbital_core/Types.h"
#include <iostream>
#include <stdexcept>

#include "EphemeridesContainer.h"
#include "earth_orientation/EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "obs/ObserverGeometry.h"

using namespace Convertors;

ObsNode::ObsNode()
    : JD{0.0},
      utcMjdDay{0},
      utcSecondsOfDay{0.0},
      hasSplitUtc{false},
      RA{0.0},
      DEC{0.0},
      MAG{0.0},
      longitude{0.0},
      cosPhi{0.0},
      sinPhi{0.0},
      LST{0.0},
      note1{' '},
      note2{' '},
      isInterstellar{false},
      magBand{' '},
      magString{},
      catalogCode{' '},
      timeAccuracyDays{1.0e-6},
      raAccuracyArcsec{1.5e-2},
      decAccuracyArcsec{1.0e-2},
      obsName{},
      vector{0.0, 0.0, 0.0},
      obsVector{0.0, 0.0, 0.0}
{
}

ObsNode::ObsNode(const orbital_core::obs::ParsedObservation& data, const ObservatoryRecord& observatory)
    : id{data.id.empty() ? orbital_core::newObservationId() : data.id},
      JD{data.jd},
      timeScale{data.timeScale},
      utcMjdDay{data.utcMjdDay},
      utcSecondsOfDay{data.utcSecondsOfDay},
      hasSplitUtc{data.hasSplitUtc},
      RA{0.0},
      DEC{0.0},
      MAG{0.0},
      longitude{observatory.longitude},
      cosPhi{observatory.cosPhi},
      sinPhi{observatory.sinPhi},
      LST{0.0},
      note1{data.programCode != ' ' ? data.programCode : (data.notes.empty() ? ' ' : data.notes.front())},
      note2{data.rovingObserver ? 'V'
            : data.isSatellite && !orbital_core::obs::techniqueOverridesGeometry(data.technique)
                ? 'S'
                : orbital_core::obs::mpcNote2(data.technique)},
      programCode{data.programCode},
      adesProgramCode{data.adesProgramCode},
      notes{data.notes},
      isInterstellar{data.isInterstellar},
      magBand{data.magBand},
      nuclearMagnitude{data.nuclearMagnitude},
      magString{data.magnitudeText},
      rmsMag{data.rmsMag},
      catalogCode{data.catalogCode},
      astCat{data.astCat},
      astrometricConvention{data.astrometricConvention},
      hasTimeAccuracy{data.hasTimeAccuracy},
      timeAccuracyDays{data.timeAccuracyDays},
      raAccuracyArcsec{data.raAccuracyArcsec},
      decAccuracyArcsec{data.decAccuracyArcsec},
      obsCode{observatory.code},
      obsName{observatory.name},
      vector{0.0, 0.0, 0.0},
      obsVector{0.0, 0.0, 0.0}
{
	if (observatory.code != data.obsCode) {
		throw std::string("Observatory record mismatch: expected '") + data.obsCode + "', got '" + observatory.code +
		    "'";
	}

	importedFit = data.importedFit;
	rovingObserver = data.rovingObserver;
	if (rovingObserver) {
		if (data.isSatellite || data.radar)
			throw std::invalid_argument("Conflicting observer geometry");
		const auto observer = orbital_core::obs::rovingObserverData(*rovingObserver, obsCode);
		longitude = observer.longitude;
		cosPhi = observer.cosPhi;
		sinPhi = observer.sinPhi;
	}
	radar = data.radar;
	if (radar) {
		if (!radar->valid())
			throw std::string("Invalid radar observation");
		note2 = 'R';
		RA = DEC = std::numeric_limits<double>::quiet_NaN();
		return;
	}

	hasObserverSigma = data.hasObserverSigma;
	rmsRaArcsec = data.rmsRaArcsec;
	rmsDecArcsec = data.rmsDecArcsec;
	rmsCorr = data.rmsCorr;
	hasRmsTime = data.hasRmsTime;
	rmsTimeSec = data.rmsTimeSec;
	hasUncTime = data.hasUncTime;
	uncTimeSec = data.uncTimeSec;

	RA = data.raRad;
	DEC = data.decRad;

	const double mjdUtc = hasSplitUtc ? static_cast<double>(utcMjdDay) + utcSecondsOfDay / 86400.0 : JD - 2400000.5;
	const double mjdTdt = mjdUtc + ttOffsetSeconds() / 86400.0;
	const double jdTdt = 2400000.5 + mjdTdt;
	double mjdUt1 = mjdUtc;
	(void)orbital_core::earthOrientation().ttToUt1Mjd(mjdTdt, mjdUt1);

	const double gast =
	    EarthRotation::Iau1976_1980::gmstRadiansMjd(mjdUt1) + EarthRotation::Iau1976_1980::equationOfEquinoxes(jdTdt);
	LST = Helpers::wrapTwoPi(gast + Convertors::degreesToRadians(longitude));

	vector = Vec3D(cos(DEC) * cos(RA), cos(DEC) * sin(RA), sin(DEC));

	MAG = stringToDouble(data.magnitudeText);

	if (data.isSatellite) {
		const double scale = (data.offsetUnits == orbital_core::obs::OffsetUnits::Au) ? 1.0 : (1.0 / 149597870.7);
		obsVector = Vec3D(-data.offsetX * scale, -data.offsetY * scale, -data.offsetZ * scale);

		isSatellite = true;
		satelliteOffsetEclAu = Convertors::ecvToEcl(
		    Vec3D(data.offsetX * scale, data.offsetY * scale, data.offsetZ * scale), Constants::EARTH_INCL);
	} else {
		obsVector =
		    Vec3D(-Constants::EARTH_RADIUS_IN_AU * cosPhi * std::cos(LST),
		          -Constants::EARTH_RADIUS_IN_AU * cosPhi * std::sin(LST), -Constants::EARTH_RADIUS_IN_AU * sinPhi);

		EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTdt, obsVector, -1);
	}

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	const auto& R = eph.getMjd(epochMjdTdb()).R;

	if (eph.error() != EphemeridesContainer::Error::NoError) {
		throw eph.errorString();
		return;
	}

	obsVector = obsVector + R;

	vector = Convertors::ecvToEcl(vector, Constants::EARTH_INCL);
	obsVector = Convertors::ecvToEcl(obsVector, Constants::EARTH_INCL);
}

double ObsNode::epochMjdTdb() const
{
	if (hasSplitUtc) {
		const double tt = static_cast<double>(utcMjdDay) + (utcSecondsOfDay + ttOffsetSeconds()) / 86400.0;
		return orbital_core::time::fromMjdTt(tt, orbital_core::time::Scale::TDB);
	}
	const double tt = orbital_core::time::mjdFromJd(JD) + ttOffsetSeconds() / 86400.0;
	return orbital_core::time::fromMjdTt(tt, orbital_core::time::Scale::TDB);
}

double ObsNode::ttOffsetSeconds() const
{
	const double mjd =
	    hasSplitUtc ? static_cast<double>(utcMjdDay) + utcSecondsOfDay / 86400.0 : orbital_core::time::mjdFromJd(JD);
	const double seconds = orbital_core::time::ttOffsetSeconds(mjd, timeScale);
	if (!std::isfinite(seconds))
		throw std::runtime_error("Observation time conversion unavailable at JD " + std::to_string(JD) +
		                         (timeScale == orbital_core::time::Scale::UT1 ? "; check earth_orientation/ET-UT.dat"
		                                                                      : "; invalid UTC epoch"));
	return seconds;
}

void ObsNode::print(std::ostream& os) const noexcept
{
	os << std::fixed << std::setprecision(6);
	os << "JD: " << JD << "\n";
	os << "obsName: " << obsName << "\n";
	os << "longitude: " << longitude << " deg\n";
	os << "MAG: " << MAG << "\n";
	os << "RA: " << RA << " rad\n";
	os << "DEC: " << DEC << " rad\n";
	os << "cosPhi: " << cosPhi << " rad\n";
	os << "sinPhi: " << sinPhi << " rad\n";
	os << "LST: " << LST << " rad, " << Convertors::radiansToDegrees(LST) << " deg, "
	   << Convertors::radiansToDegrees(LST) / 15 << " hours\n";
	os.unsetf(std::ios::fixed);
	os << "\nvector (L): " << vector << "\n";
	os << "obsVector (R): " << obsVector << "\n";
}
