#include "AstrometryBias.h"
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstring>
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"
#include "io/DataSourceFactory.h"
#include "logger/Logger.h"

#include <array>
#include <mutex>
#include <unordered_map>

namespace {
	constexpr int kNpix = 12 * AstrometryBias::NSIDE * AstrometryBias::NSIDE;
	constexpr int kValsPerCat = 4;

	constexpr std::array<char, 26> kFarnocchia2018CatalogOrder = {'a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l',
	                                                              'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v',
	                                                              'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'Y'};

	constexpr std::array<char, 19> kFcct14CatalogOrder = {'a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm',
	                                                      'o', 'p', 'q', 'r', 'u', 'v', 'w', 'L', 'N'};

	constexpr double kDaysPerJulianYear = 365.25;

	constexpr double kBiasReferenceJd = 2'451'544.5;
	constexpr double kMasPerArcsec = 1000.0;
	constexpr double kCosDecEpsilon = 1e-12;
	constexpr double kHealpixEquatorialZLimit = 2.0 / 3.0;

	double computeDeltaYears(double jd) { return (jd - kBiasReferenceJd) / kDaysPerJulianYear; }

	double safeCosDec(double dec)
	{
		double c = std::cos(dec);
		if (std::fabs(c) < kCosDecEpsilon)
			return (c < 0.0) ? -kCosDecEpsilon : kCosDecEpsilon;

		return c;
	}

	Vec3D toUnitVector(double ra, double dec)
	{
		const double cosDec = std::cos(dec);
		return {cosDec * std::cos(ra), cosDec * std::sin(ra), std::sin(dec)};
	}

	AstrometryBias& getBias(const std::filesystem::path& path)
	{
		static std::mutex mutex;
		static std::unordered_map<std::string, std::unique_ptr<AstrometryBias>> cache;

		const std::string key = path.lexically_normal().string();
		std::lock_guard lock(mutex);
		const auto found = cache.find(key);
		if (found != cache.end())
			return *found->second;

		auto bias = std::make_unique<AstrometryBias>();
		try {
			bias->load(path);
		} catch (const std::exception& e) {
			LOG_INFO << "AstrometryBias load failed: " << e.what();
		}

		auto* result = bias.get();
		cache.emplace(key, std::move(bias));
		return *result;
	}

}  // namespace

class AstrometryBias::Impl {
   public:
	std::vector<double> data;
	std::array<int, 256> catMap{};
	int rowColumns{};

	Impl() { catMap.fill(-1); }

	[[nodiscard]] int catalogIndex(char code) const noexcept { return catMap[static_cast<unsigned char>(code)]; }

	[[nodiscard]] std::size_t cellBase(int pix, int cat) const noexcept
	{
		return static_cast<std::size_t>(pix) * static_cast<std::size_t>(rowColumns) +
		       static_cast<std::size_t>(cat) * kValsPerCat;
	}
};

AstrometryBias::AstrometryBias() : p_impl(std::make_unique<Impl>()) {}

AstrometryBias::~AstrometryBias() = default;
AstrometryBias::AstrometryBias(AstrometryBias&&) noexcept = default;
AstrometryBias& AstrometryBias::operator=(AstrometryBias&&) noexcept = default;

void AstrometryBias::load(const std::filesystem::path& path)
{
	auto source = DataSourceFactory::open(path.string());
	if (!source || !source->isOpen()) {
		throw std::runtime_error("AstrometryBias: cannot open " + path.string());
	}

	const auto bytes = source->bytes();

	std::vector<char> buf(bytes.size());
	if (!bytes.empty()) {
		std::memcpy(buf.data(), bytes.data(), bytes.size());
	}

	const char* p = buf.data();
	const char* end = p + buf.size();

	auto isWs = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

	auto skipToData = [&](const char*& cursor) {
		while (cursor < end) {
			while (cursor < end && (isWs(*cursor) || *cursor == '\n'))
				++cursor;

			if (cursor >= end || *cursor != '!')
				return;

			while (cursor < end && *cursor != '\n')
				++cursor;
		}
	};

	skipToData(p);
	if (p >= end)
		throw std::runtime_error("AstrometryBias: empty table " + path.string());

	const char* firstLineEnd = p;
	while (firstLineEnd < end && *firstLineEnd != '\n')
		++firstLineEnd;

	int firstRowColumns = 0;
	for (const char* q = p; q < firstLineEnd;) {
		while (q < firstLineEnd && isWs(*q))
			++q;
		if (q >= firstLineEnd)
			break;
		double value = 0.0;
		const auto parsed = std::from_chars(q, firstLineEnd, value);
		if (parsed.ec != std::errc{})
			throw std::runtime_error("AstrometryBias: invalid first data row in " + path.string());

		++firstRowColumns;
		q = parsed.ptr;
	}

	p_impl->catMap.fill(-1);
	if (firstRowColumns == static_cast<int>(kFarnocchia2018CatalogOrder.size()) * kValsPerCat) {
		for (std::size_t i = 0; i < kFarnocchia2018CatalogOrder.size(); ++i) {
			p_impl->catMap[static_cast<unsigned char>(kFarnocchia2018CatalogOrder[i])] = static_cast<int>(i);
		}
	}

	else if (firstRowColumns == static_cast<int>(kFcct14CatalogOrder.size()) * kValsPerCat) {
		for (std::size_t i = 0; i < kFcct14CatalogOrder.size(); ++i) {
			p_impl->catMap[static_cast<unsigned char>(kFcct14CatalogOrder[i])] = static_cast<int>(i);
		}

		p_impl->catMap[static_cast<unsigned char>('g')] = -1;
		p_impl->catMap[static_cast<unsigned char>('l')] = -1;
	}

	else {
		throw std::runtime_error("AstrometryBias: unsupported column count " + std::to_string(firstRowColumns) +
		                         " in " + path.string());
	}

	p_impl->rowColumns = firstRowColumns;
	p_impl->data.assign(static_cast<std::size_t>(kNpix) * static_cast<std::size_t>(p_impl->rowColumns), 0.0);

	int pix = 0;

	while (p < end) {
		skipToData(p);
		if (p >= end)
			break;

		if (pix >= kNpix)
			throw std::runtime_error("Too many rows in astrometric bias table");

		const char* lineEnd = p;
		while (lineEnd < end && *lineEnd != '\n')
			++lineEnd;

		const std::size_t base = static_cast<std::size_t>(pix) * static_cast<std::size_t>(p_impl->rowColumns);

		for (int k = 0; k < p_impl->rowColumns; ++k) {

			while (p < lineEnd && isWs(*p))
				++p;

			if (p >= lineEnd)
				throw std::runtime_error("Truncated row at pix=" + std::to_string(pix));

			double v = 0.0;

			auto res = std::from_chars(p, lineEnd, v);
			if (res.ec != std::errc{})
				throw std::runtime_error("Parse error at pix=" + std::to_string(pix) + " col=" + std::to_string(k));

			p_impl->data[base + k] = v;
			p = res.ptr;
		}

		++pix;
		p = lineEnd;
	}

	if (pix != kNpix)
		throw std::runtime_error("Wrong row count in astrometric bias table");

	LOG_INFO << "AstrometryBias loaded: " << pix << " rows, " << p_impl->rowColumns / kValsPerCat << " catalogs from "
	         << path.string();
}

bool AstrometryBias::loaded() const noexcept { return !p_impl->data.empty(); }

std::size_t AstrometryBias::rowCount() const noexcept
{
	return loaded() && p_impl->rowColumns > 0 ? p_impl->data.size() / p_impl->rowColumns : 0;
}

std::size_t AstrometryBias::catalogCount() const noexcept
{
	return loaded() ? static_cast<std::size_t>(p_impl->rowColumns / kValsPerCat) : 0;
}

AstrometryBias::Correction AstrometryBias::correction(double raRad, double decRad, char catalogCode, double jd) const
{
	if (!loaded())
		return {};

	const int cat = p_impl->catalogIndex(catalogCode);
	if (cat < 0)
		return {};

	const double theta = Constants::PI * 0.5 - decRad;
	const int pix = ang2pixRing(theta, raRad, NSIDE);

	if (pix < 0 || pix >= kNpix)
		return {};

	const std::size_t base = p_impl->cellBase(pix, cat);

	const double dRA_cos = p_impl->data[base + 0];
	const double dDEC = p_impl->data[base + 1];
	const double pmRA_cos = p_impl->data[base + 2];
	const double pmDEC = p_impl->data[base + 3];

	const double dt = computeDeltaYears(jd);
	const double cosDec = safeCosDec(decRad);

	return {(dRA_cos + dt * pmRA_cos / kMasPerArcsec) / cosDec, dDEC + dt * pmDEC / kMasPerArcsec};
}

int AstrometryBias::ang2pixRing(double theta, double phi, int nside) noexcept
{
	const double z = std::cos(theta);
	const double za = std::fabs(z);

	double phi_mod = std::fmod(phi, 2.0 * Constants::PI);
	if (phi_mod < 0.0)
		phi_mod += 2.0 * Constants::PI;

	const double tt = phi_mod / (0.5 * Constants::PI);

	int pix = -1;

	if (za <= kHealpixEquatorialZLimit) {
		const double temp1 = nside * (0.5 + tt);
		const double temp2 = nside * 0.75 * z;

		const int jp = static_cast<int>(std::floor(temp1 - temp2));
		const int jm = static_cast<int>(std::floor(temp1 + temp2));

		const int ir = nside + 1 + jp - jm;
		const int kshift = 1 - (ir & 1);

		int ip = (jp + jm - nside + kshift + 1) / 2;
		ip = ((ip % (4 * nside)) + 4 * nside) % (4 * nside);

		pix = 2 * nside * (nside - 1) + (ir - 1) * 4 * nside + ip;
	}

	else {
		const double tp = tt - std::floor(tt);
		const double tmp = nside * std::sqrt(3.0 * (1.0 - za));

		const int jp = static_cast<int>(std::floor(tp * tmp));
		const int jm = static_cast<int>(std::floor((1.0 - tp) * tmp));

		const int ir = jp + jm + 1;

		int ip = static_cast<int>(std::floor(tt * ir));
		ip = ((ip % (4 * ir)) + 4 * ir) % (4 * ir);

		if (z > 0.0)
			pix = 2 * ir * (ir - 1) + ip;

		else
			pix = 12 * nside * nside - 2 * ir * (ir + 1) + ip;
	}

	return pix;
}

std::filesystem::path Calculations::astrometryBiasPathForScheme(const std::filesystem::path& runtimeDirectory,
                                                                AstrometryBiasScheme scheme)
{
	const bool efcc18 = scheme == AstrometryBiasScheme::Efcc18;
	const auto current = runtimeDirectory / "debias" / (efcc18 ? "efcc18.dat" : "fcct14.dat");
	const auto legacy = efcc18 ? runtimeDirectory / "debias_2018" / "bias.dat"
	                           : runtimeDirectory / "debias_fcct14" / "xcat.bias_fcct14";
	std::error_code ec;
	if (std::filesystem::is_regular_file(current, ec))
		return current;
	ec.clear();
	if (std::filesystem::is_regular_file(legacy, ec))
		return legacy;
	return current;
}

std::filesystem::path Calculations::astrometryBiasPathForProfile(const std::filesystem::path& runtimeDirectory,
                                                                 orbital_core::AstrometricWeightProfile /*profile*/)
{
	return astrometryBiasPathForScheme(runtimeDirectory, AstrometryBiasScheme::Efcc18);
}

std::size_t Calculations::debiasObservations(std::vector<ObsNode>& obs, const std::filesystem::path& biasPath)
{
	auto& bias = getBias(biasPath);
	if (!bias.loaded())
		return 0;

	constexpr double kArcsecToRad = Constants::PI / (180.0 * 3600.0);
	std::size_t corrected = 0;

	for (auto& o : obs) {
		if (o.radar)
			continue;
		const double rawRa = o.submittedRaRad();
		const double rawDec = o.submittedDecRad();
		const auto c = bias.correction(rawRa, rawDec, o.catalogCode, o.JD);

		if (c.dRA_arcsec == 0.0 && c.dDEC_arcsec == 0.0)
			continue;

		if (!o.submittedCoordinates)
			o.submittedCoordinates = ObsNode::Coordinates{rawRa, rawDec, o.vector};

		o.RA = rawRa - c.dRA_arcsec * kArcsecToRad;
		o.DEC = rawDec - c.dDEC_arcsec * kArcsecToRad;
		o.biasRaArcsec = c.dRA_arcsec;
		o.biasDecArcsec = c.dDEC_arcsec;

		const Vec3D vEq = toUnitVector(o.RA, o.DEC);
		o.vector = Convertors::ecvToEcl(vEq, Constants::EARTH_INCL);

		++corrected;
	}

	return corrected;
}
