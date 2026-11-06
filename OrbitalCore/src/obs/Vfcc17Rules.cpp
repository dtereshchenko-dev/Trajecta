#include "obs/Vfcc17Rules.h"
#include <string>
#include <string_view>
#include <vector>

#include <array>
#include <charconv>
#include <chrono>
#include <utility>

namespace orbital_core::obs {
	namespace {

		std::string_view trim(std::string_view s) noexcept
		{
			const auto first = s.find_first_not_of(" \t\r");
			return first == s.npos ? std::string_view{} : s.substr(first, s.find_last_not_of(" \t\r") - first + 1);
		}

		template <typename T>
		bool number(std::string_view s, T& value) noexcept
		{
			s = trim(s);
			if (s.empty())
				return false;
			const auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
			return ec == std::errc{} && end == s.data() + s.size();
		}

		bool date(std::string_view s, double& jd) noexcept
		{
			s = trim(s);
			if (s.empty())
				return true;
			int y = 0, m = 0, d = 0;
			if (s.size() != 10 || s[4] != '-' || s[7] != '-' || !number(s.substr(0, 4), y) ||
			    !number(s.substr(5, 2), m) || !number(s.substr(8, 2), d))
				return false;
			const std::chrono::year_month_day civil{std::chrono::year{y}, std::chrono::month{static_cast<unsigned>(m)},
			                                        std::chrono::day{static_cast<unsigned>(d)}};
			if (!civil.ok())
				return false;
			jd = 2440587.5 + static_cast<double>(std::chrono::sys_days{civil}.time_since_epoch().count());
			return true;
		}

	}  // namespace

	bool Vfcc17Rules::load(std::string_view text, std::string& error)
	{
		error.clear();
		m_rules.clear();
		std::vector<Rule> parsed;
		std::size_t lineNumber = 0;
		if (text.starts_with("\xEF\xBB\xBF"))
			text.remove_prefix(3);
		while (!text.empty()) {
			const auto end = text.find('\n');
			auto line = text.substr(0, end);
			text = end == text.npos ? std::string_view{} : text.substr(end + 1);
			++lineNumber;

			const auto content = trim(line);
			if (content.empty() || content.front() == '!' || content.front() == '#')
				continue;

			const auto fail = [&] {
				error = "vfcc17.rules: malformed rule at line " + std::to_string(lineNumber);
				return false;
			};

			if (line.size() < 71 || line.substr(4, 2) != "t=" || line.substr(12, 2) != "c=" ||
			    line.substr(25, 2) != "p=" || line[31] != '>' || line[44] != '<' || line[57] != '@' || line[64] != ',')
				return fail();

			Rule r;
			r.station = trim(line.substr(0, 3));
			r.techniques = trim(line.substr(6, 5));
			r.catalogs = trim(line.substr(14, 10));
			r.programs = trim(line.substr(27, 3));
			if (r.station.size() != 3 || r.techniques.empty() || r.catalogs.empty() ||
			    !date(line.substr(33, 10), r.startJd) || !date(line.substr(46, 10), r.endJd) || r.startJd >= r.endJd ||
			    !number(line.substr(59, 5), r.sigmaRa) || !number(line.substr(66, 5), r.sigmaDec) ||
			    !(r.sigmaRa > 0.0 && r.sigmaRa < 1.0e6) || !(r.sigmaDec > 0.0 && r.sigmaDec < 1.0e6))
				return fail();
			parsed.push_back(std::move(r));
		}

		if (parsed.empty()) {
			error = "vfcc17.rules: no rules";
			return false;
		}

		m_rules = std::move(parsed);
		return true;
	}

	size_t Vfcc17Rules::size() const noexcept { return m_rules.size(); }

	ObservationUncertainty Vfcc17Rules::lookup(std::string_view station, char technique, char catalog, char program,
	                                           double jdUtc) const noexcept
	{
		ObservationUncertainty out;
		if (technique == 'x' || technique == 'X') {
			out.sigmaRaRad = out.sigmaDecRad = 999.0 * Constants::ARCSEC_RAD;
			out.resolved = true;
			return out;
		}

		constexpr std::array<std::array<bool, 3>, 5> passes{
		    {{true, true, true}, {true, true, false}, {true, false, true}, {true, false, false}, {false, false, true}}};

		for (const auto scope : {station, std::string_view{"ALL"}}) {
			for (const auto& pass : passes) {
				for (const auto& r : m_rules) {
					if (r.station != scope || !(jdUtc > r.startJd && jdUtc <= r.endJd))
						continue;
					if (r.techniques.find(pass[0] ? technique : '*') == std::string::npos ||
					    r.catalogs.find(pass[1] ? catalog : '*') == std::string::npos)
						continue;
					if (pass[2] ? (program == ' ' ? !r.programs.empty() : r.programs.find(program) == std::string::npos)
					            : !r.programs.empty())
						continue;
					out.sigmaRaRad = r.sigmaRa * Constants::ARCSEC_RAD;
					out.sigmaDecRad = r.sigmaDec * Constants::ARCSEC_RAD;
					out.resolved = true;
					return out;
				}
			}
		}
		return out;
	}

}  // namespace orbital_core::obs
