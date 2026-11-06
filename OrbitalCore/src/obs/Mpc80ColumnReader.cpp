#include "obs/Mpc80ColumnReader.h"
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "MPECParser.h"
#include "orbital_core/Types.h"

#include <cmath>
#include <sstream>

namespace orbital_core::obs {
	namespace {

		[[nodiscard]] Technique techniqueFromNote2(char note2) noexcept
		{
			switch (note2) {
				case 'C':
					return Technique::Ccd;
				case 'c':
					return Technique::CorrectedCcd;
				case 'B':
					return Technique::Cmos;
				case ' ':
				case 'P':
					return Technique::Photographic;
				case 'A':
					return Technique::PhotographicB1950;
				case 'H':
					return Technique::Photometer;
				case 'X':
					return Technique::ReplacedCcd;
				case 'x':
					return Technique::ReplacedPhotographic;
				case 'n':
					return Technique::Video;
				case 'e':
					return Technique::Encoder;
				case 'M':
					return Technique::Micrometer;
				case 'T':
					return Technique::Meridian;
				case 'E':
					return Technique::Occultation;
				case 'N':
					return Technique::NormalPlace;
				default:
					return Technique::Unknown;
			}
		}

		[[nodiscard]] bool sexagesimalToRadians(const std::string& raText, const std::string& decText, double& raRad,
		                                        double& decRad) noexcept
		{
			if (raText.size() < 7 || decText.size() < 8)
				return false;
			try {
				const int raH = Convertors::stringToInt(raText.substr(0, 2));
				const int raM = Convertors::stringToInt(raText.substr(3, 2));
				const double raS = Convertors::stringToDouble(raText.substr(6));
				raRad = Convertors::raToRadians(raH, raM, raS);

				const bool negative = decText.front() == '-';
				const int decD = std::abs(Convertors::stringToInt(decText.substr(0, 3)));
				const int decM = Convertors::stringToInt(decText.substr(4, 2));
				const double decS = Convertors::stringToDouble(decText.substr(7));
				decRad = Convertors::decToRadians(decD, decM, decS);
				if (negative)
					decRad = -decRad;
				return true;
			} catch (...) {
				return false;
			}
		}

	}  // namespace

	const char* Mpc80ColumnReader::formatName() const noexcept { return "MPC 80-column"; }

	bool Mpc80ColumnReader::canRead(std::string_view text) const noexcept
	{
		const std::size_t scanEnd = std::min<std::size_t>(text.size(), 65536);
		std::size_t begin = 0;
		while (begin < scanEnd) {
			std::size_t end = text.find('\n', begin);
			if (end == std::string_view::npos || end > scanEnd)
				end = scanEnd;

			const std::string_view line = text.substr(begin, end - begin);
			begin = end + 1;
			if (line.empty())
				continue;

			if (line.size() >= 80 && (line[44] == '+' || line[44] == '-' || line[14] == 'R' || line[14] == 'Q'))
				return true;

			std::size_t fields = 0;
			bool inField = false;
			for (const char c : line) {
				const bool space = (c == ' ' || c == '\t' || c == '\r');
				if (!space && !inField) {
					++fields;
					inField = true;
				}

				else if (space) {
					inField = false;
				}
			}

			if (fields >= 10)
				return true;
		}

		return false;
	}

	std::vector<ParsedObservation> Mpc80ColumnReader::read(const std::string& text, std::string& error) const
	{
		error.clear();

		std::vector<std::string> lines;
		std::istringstream stream(text);
		std::string line;
		while (std::getline(stream, line))
			lines.push_back(line);

		MPECParser parser(lines);
		const std::vector<MPECInitData> raw = parser.parse();
		if (parser.error() != MPECParser::Error::NoError) {
			error = parser.errorString();
			return {};
		}

		std::vector<ParsedObservation> out;
		out.reserve(raw.size());
		for (const MPECInitData& d : raw) {
			ParsedObservation o;
			o.jd = d.JD;
			o.timeScale = time::legacyObservationScale(o.jd);
			o.utcMjdDay = d.utcMjdDay;
			o.utcSecondsOfDay = d.utcSecondsOfDay;
			o.hasSplitUtc = d.hasSplitUtc;

			o.radar = d.radar;
			if (!o.radar && !sexagesimalToRadians(d.ra, d.dec, o.raRad, o.decRad)) {
				error = "MPC 80-column: malformed RA/DEC field";
				return {};
			}

			o.objectName = d.objectName;
			o.classificationDesignation = d.classificationDesignation;
			o.obsCode = d.code;
			o.technique = techniqueFromNote2(d.note2);
			o.programCode = d.note1;
			o.isInterstellar = d.isInterstellar;

			if (!d.catalog.empty())
				o.catalogCode = d.catalog[0];
			else if (d.magBand.empty() && !d.filter.empty())
				o.catalogCode = d.filter[0];

			o.magnitudeText = d.magnitude;
			o.magBand = d.magBand.empty() ? ' ' : d.magBand[0];
			if (o.magBand == 'N' || o.magBand == 'T')
				o.nuclearMagnitude = o.magBand == 'N';

			o.timeAccuracyDays = d.timeAccuracyDays;
			o.hasTimeAccuracy = d.timeAccuracyDays > 0.0;
			o.raAccuracyArcsec = d.raAccuracyArcsec * std::max(std::abs(std::cos(o.decRad)), 1.0e-12);
			o.decAccuracyArcsec = d.decAccuracyArcsec;

			o.isSatellite = d.isSatellite;
			o.rovingObserver = d.rovingObserver;
			if (d.isSatellite) {
				o.offsetUnits = d.satOffsetUnits == '1' ? OffsetUnits::Km : OffsetUnits::Au;
				o.offsetX = d.satOffsetX;
				o.offsetY = d.satOffsetY;
				o.offsetZ = d.satOffsetZ;
			}

			o.sourceFormat = "MPC 80-column";
			out.push_back(std::move(o));
		}

		if (out.empty())
			error = "MPC 80-column: no observations parsed";
		return out;
	}

}  // namespace orbital_core::obs
