#ifndef MPECPARSER_H
#define MPECPARSER_H
#include <cstdint>
#include <string>
#include <vector>
#include "MPECInitData.h"
#include "orbital_core/OrbitalCoreExport.h"

class ORBITAL_CORE_API MPECParser {
   public:
	enum class Error : int {
		NoError = 0,
		InvalidData,
		InsufficientData,
		InvalidRAData,
		InvalidDecData,
		InvalidMagnitude,
		InvalidObjectName,  // Новая ошибка: Некорректное имя объекта
		InvalidDate,        // Новая ошибка: Ошибка преобразования даты
		InvalidCode,        // Новая ошибка: Ошибка кода обсерватории
		InvalidFormat       // Новая ошибка: Ошибка формата данных
	};

	MPECParser(const std::vector<std::string>& lines);
	~MPECParser() = default;

	std::vector<MPECInitData> parse() noexcept;

	Error error() const noexcept;
	std::string errorString() const noexcept;

   private:
	MPECParser() = delete;
	MPECParser(const MPECParser&) = delete;
	MPECParser(MPECParser&&) = delete;
	MPECParser& operator=(const MPECParser&) = delete;
	MPECParser& operator=(MPECParser&&) = delete;

	MPECInitData parseLine(const std::string& line);

	// Strict MPC 80-column optical-observation format (official MPC submission).
	// Returns a fully-populated MPECInitData or sets m_error.
	MPECInitData parseLine80Col(const std::string& line);

	// Legacy MPC radar dates encode whole UTC seconds as rounded decimal days.
	// Recover those seconds and normalize the split MJD across midnight; this
	// decoding rule does not apply to optical MPC, ADES or RWO timestamps.
	bool parseRadarPair(const std::string& first, const std::string& second, MPECInitData& data);

	// Legacy whitespace-tokenised parser (kept for user-supplied short lines
	// where the 80-column layout is not enforced).
	MPECInitData parseLineTokenised(const std::string& line);

	// True when the line conforms to strict MPC 80-column layout:
	// at least 80 chars AND a Dec sign at column 45 (0-indexed position 44).
	static bool isMPC80Col(const std::string& line) noexcept;

	// True when the line is a satellite-observation continuation line —
	// 's' (lowercase) at column 15 (0-indexed 14).  These come AFTER an
	// ordinary obs line whose note-2 is 'S' (uppercase) and carry the
	// spacecraft → geocenter offset in J2000 equatorial X/Y/Z.
	static bool isSatelliteContinuation(const std::string& line) noexcept;

	// True when the line is a roving-observer continuation line —
	// 'v' (lowercase) at column 15 (0-indexed 14).  These come AFTER an
	// ordinary obs line whose note-2 is 'V' (uppercase) and carry the
	// observer's East longitude / geodetic latitude / altitude.  MPC
	// uses obs code "247" on the parent line as a sentinel.
	static bool isRovingContinuation(const std::string& line) noexcept;

	// Parse the satellite continuation line and write offsets into `into`.
	// Returns false on malformed input.  Units are inferred from col 33:
	//   '1' = AU, '2' = km.  Sign + magnitude tokens are read by
	//   whitespace tokenisation of the substring after col 32 to be
	//   tolerant of MPC's irregular spacing in real-world feeds.
	bool parseSatelliteContinuation(const std::string& line, MPECInitData& into) noexcept;

	// Parse the roving-observer continuation line and write the equivalent
	// geocenter→observer offset into `into` (km, J2000 equatorial).  The
	// raw input is (E.lon, lat, alt) geodetic; we convert geodetic→ECEF
	// (WGS84 ellipsoid) then ECEF→ECI via GMST(JD).  Stored using the
	// same satOffsetX/Y/Z + satOffsetUnits='1' (km) fields as satellite,
	// and isSatellite=true marks the parent — downstream geometry is
	// identical (observer position read from continuation, not from
	// observatory-code geodetic table).
	bool parseRovingContinuation(const std::string& line, MPECInitData& into) noexcept;

   private:
	bool isDouble(const std::string& s) noexcept;
	bool isNumber(const std::string& s) noexcept;

   private:
	std::vector<std::string> m_lines;
	Error m_error = Error::NoError;

	constexpr static const char RA_DELIMITER = ' ';
	constexpr static const char DEC_DELIMITER = ' ';
	constexpr static const char DEC_SIGN = '+';
	constexpr static const char DEC_SIGN_NEGATIVE = '-';
	constexpr static uint8_t MIN_TOKENS_COUNT = 10;
};

#endif  // MPECPARSER_H
