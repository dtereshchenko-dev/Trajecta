#include "orbital_core/Residuals.h"
#include <filesystem>
#include <string>
#include <vector>
#include "orbital_core/Types.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

void Residuals::add(const Residual& r) { m_residuals.push_back(r); }

void Residuals::clear() { m_residuals.clear(); }

double Residuals::rms() const
{
	if (m_residuals.empty())
		return 0.0;

	double sumSq = 0.0;
	std::size_t opticalCount = 0;
	for (const auto& r : m_residuals) {
		if (r.radar)
			continue;
		sumSq += r.dRA_arcsec * r.dRA_arcsec + r.dDEC_arcsec * r.dDEC_arcsec;
		++opticalCount;
	}

	return opticalCount ? std::sqrt(sumSq / (2.0 * opticalCount)) : 0.0;
}

size_t Residuals::count() const noexcept { return m_residuals.size(); }

const std::vector<Residual>& Residuals::data() const noexcept { return m_residuals; }

void Residuals::printHeader(std::ostream& os) const
{
	os << std::left << std::setw(14) << "Time" << std::setw(6) << "Code" << std::setw(12) << " RA" << std::setw(14)
	   << "DEC" << std::right << std::setw(7) << "RMSRa" << std::setw(8) << "RMSDec" << std::setw(8) << "Total"
	   << "\n";
}

void Residuals::print() const { print(std::cout); }

void Residuals::print(std::ostream& os) const
{
	printHeader(os);

	auto fmtResid = [](double value) -> std::string {
		char sign = (value >= 0.0) ? '+' : '-';
		double abs = std::fabs(value);

		std::ostringstream oss;
		oss << std::fixed << std::setprecision(2) << abs;
		return oss.str() + sign;
	};

	os << std::fixed;
	for (const auto& r : m_residuals) {
		if (r.radar) {
			os << Convertors::jdToPackedDate(r.jd) << " Radar " << r.radar->transmitter << "->" << r.obsCode;
			if (r.radar->delaySeconds)
				os << " delay O-C=" << std::setprecision(9) << r.radar->delaySeconds->residual << " s";
			if (r.radar->dopplerHz)
				os << " Doppler O-C=" << std::setprecision(6) << r.radar->dopplerHz->residual << " Hz";
			os << '\n';
			continue;
		}
		int raH, raM;
		double raS;
		Convertors::radiansToHMS(r.raObs, raH, raM, raS);
		int decD, decM;
		double decS;
		Convertors::radiansToDMS(r.decObs, decD, decM, decS);

		std::ostringstream raStr;
		raStr << std::fixed << std::setw(2) << raH << ":" << std::setfill('0') << std::setw(2) << raM << ":"
		      << std::setprecision(2) << std::setw(5) << raS << std::setfill(' ');

		std::ostringstream decStr;
		decStr << std::fixed << std::showpos << decD << std::noshowpos << " " << std::setfill('0') << std::setw(2)
		       << std::abs(decM) << " " << std::setprecision(1) << std::setw(4) << std::abs(decS) << std::setfill(' ');

		os << std::left << std::setw(14) << Convertors::jdToPackedDate(r.jd) << std::setw(6) << r.obsCode
		   << std::setw(12) << raStr.str() << std::setw(14) << decStr.str() << std::right << std::setw(7)
		   << fmtResid(r.dRA_arcsec) << std::setw(7) << fmtResid(r.dDEC_arcsec) << std::setw(8) << std::fixed
		   << std::setprecision(2) << r.total_arcsec << "\n";
	}

	os << "RMS: " << std::setprecision(2) << rms() << "\" (" << m_residuals.size() << " observations)\n";
}

bool Residuals::save(const std::filesystem::path& filePath) const
{
	std::ofstream ofs(filePath);
	if (!ofs.is_open())
		return false;

	print(ofs);
	return true;
}
