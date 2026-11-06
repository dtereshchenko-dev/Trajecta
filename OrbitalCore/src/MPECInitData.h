#pragma once
#include <string>
#include <ostream>
#include <iostream>
#include "orbital_core/obs/RadarObservation.h"
#include "orbital_core/obs/RovingObserver.h"

struct MPECInitData {
	MPECInitData() = default;
	~MPECInitData() = default;

	MPECInitData(const MPECInitData& other) = default;
	MPECInitData(MPECInitData&& other) = default;
	MPECInitData& operator=(const MPECInitData& other) = default;
	MPECInitData& operator=(MPECInitData&& other) = default;

	void print(std::ostream& os = std::cout) const noexcept
	{
		os << "Object name: " << objectName << "\n";
		os << "RA: " << ra << "\n";
		os << "Dec: " << dec << "\n";
		os << "Magnitude: " << magnitude << "\n";
		os << "Filter: " << filter << "\n";
		os << "Mag band: " << magBand << "\n";
		os << "Catalog: " << catalog << "\n";
		os << "Code: " << code << "\n";
	}

	bool operator==(const MPECInitData& other) const noexcept
	{
		return JD == other.JD && objectName == other.objectName && ra == other.ra && dec == other.dec &&
		       magnitude == other.magnitude && code == other.code && radar == other.radar;
	}

	double JD;
	std::optional<orbital_core::obs::RadarObservation> radar;
	// MPC time retained in the same split representation used by reference.
	// Keeping the integer MJD separate from seconds avoids losing tens of
	// microseconds when the fractional day is added to a full Julian date.
	int utcMjdDay = 0;
	double utcSecondsOfDay = 0.0;
	bool hasSplitUtc = false;
	std::string objectName;
	std::string classificationDesignation;
	std::string ra;
	std::string dec;
	std::string magnitude;
	std::string magBand;
	std::string catalog;
	std::string filter;
	std::string code;
	char note1 = ' ';
	char note2 = ' ';
	bool isInterstellar = false;
	double timeAccuracyDays = 1.0e-6;
	double raAccuracyArcsec = 1.5e-2;
	double decAccuracyArcsec = 1.0e-2;

	// ---- Satellite-observation continuation data ----
	// Set to true when the parent obs line has note-2 = 'S' AND the following
	// line was a valid 's' continuation that we successfully parsed.
	// Offsets are spacecraft → geocenter in J2000 equatorial frame; units
	// depend on satOffsetUnits ('1' = km, '2' = AU per MPC convention).
	// ObsNode constructor converts to ecliptic AU.
	bool isSatellite = false;
	std::optional<orbital_core::obs::RovingObserver> rovingObserver;
	char satOffsetUnits = ' ';  // '1' = km, '2' = AU
	double satOffsetX = 0.0;
	double satOffsetY = 0.0;
	double satOffsetZ = 0.0;
};
