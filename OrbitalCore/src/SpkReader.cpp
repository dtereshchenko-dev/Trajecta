
#include "orbital_core/SpkReader.h"
#include "orbital_core/RuntimeData.h"
#include <memory>
#include <string>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include "orbital_core/Types.h"
#include "io/DataSourceFactory.h"
#include "io/IDataSource.h"
#include "ephemerides/SplitTime.h"
#include "logger/Logger.h"

#include <array>
#include <filesystem>
#include <unordered_map>
#include <vector>

namespace {
	constexpr int SPK_TYPE_CHEBYSHEV_POSITION = 2;
	constexpr int NAIF_SSB_ID = 0;
	constexpr int NAIF_SUN_ID = 10;
	constexpr size_t DAF_RECORD_BYTES = 1024;
	constexpr size_t SUMMARY_DOUBLES = 5;
	constexpr size_t SUMMARY_BYTES = SUMMARY_DOUBLES * sizeof(double);
	constexpr double AU_KM = 149597870.7;

	const double kCosObliquity = std::cos(Constants::EARTH_INCL);
	const double kSinObliquity = std::sin(Constants::EARTH_INCL);

	inline Vec3D rotateEqToEcl(double x, double y, double z)
	{
		return {x, y * kCosObliquity + z * kSinObliquity, z * kCosObliquity - y * kSinObliquity};
	}

	inline void swap8(void* p) noexcept
	{
		auto* b = static_cast<unsigned char*>(p);
		std::swap(b[0], b[7]);
		std::swap(b[1], b[6]);
		std::swap(b[2], b[5]);
		std::swap(b[3], b[4]);
	}

	inline void swap4(void* p) noexcept
	{
		auto* b = static_cast<unsigned char*>(p);
		std::swap(b[0], b[3]);
		std::swap(b[1], b[2]);
	}

	inline std::int32_t doubleLowInt(double v, bool byteSwap) noexcept
	{
		std::int32_t lo;
		std::memcpy(&lo, &v, 4);
		if (byteSwap)
			swap4(&lo);

		return lo;
	}

	inline std::int32_t doubleHighInt(double v, bool byteSwap) noexcept
	{
		std::int32_t hi;
		std::memcpy(&hi, reinterpret_cast<const char*>(&v) + 4, 4);
		if (byteSwap)
			swap4(&hi);

		return hi;
	}

	inline bool readDoublesAt(raw_file_data view, std::size_t offset, std::size_t n, double* out,
	                          bool byteSwap) noexcept
	{
		const std::size_t needed = n * sizeof(double);
		if (offset + needed > view.size())
			return false;

		std::memcpy(out, view.data() + offset, needed);
		if (byteSwap) {
			for (std::size_t i = 0; i < n; ++i)
				swap8(out + i);
		}

		return true;
	}

	inline bool readInt32At(raw_file_data view, std::size_t offset, std::int32_t& out, bool byteSwap) noexcept
	{
		if (offset + sizeof(out) > view.size())
			return false;
		std::memcpy(&out, view.data() + offset, sizeof(out));

		if (byteSwap)
			swap4(&out);
		return true;
	}
}  // namespace

struct SpkReader::Impl {
	std::unique_ptr<IDataSource> source;
	std::string sourcePath;
	bool available = false;
	bool byteSwap = false;

	[[nodiscard]] raw_file_data view() const noexcept
	{
		if (!source || !source->isOpen())
			return {};

		return source->bytes();
	}

	struct Segment {
		double etStart = 0.0;
		double etEnd = 0.0;
		int target = 0;
		int center = 0;
		int frame = 0;
		int dataType = 0;
		int startAddr = 0;
		int endAddr = 0;

		double initEt = 0.0;
		double intLen = 0.0;
		int rsize = 0;
		int nRecords = 0;
		int nc = 0;
		bool tailRead = false;
	};

	std::unordered_map<int, std::vector<Segment>> segByTarget;

	const Segment* pickSegment(int naifId, orbital_core::detail::SplitTime et) const
	{
		auto it = segByTarget.find(naifId);
		if (it == segByTarget.end())
			return nullptr;
		for (const auto& seg : it->second) {
			if (et.atOrAfter(seg.etStart) && et.atOrBefore(seg.etEnd))
				return &seg;
		}
		return nullptr;
	}

	bool open(const std::string& path)
	{
		sourcePath = path;
		source.reset();
		available = false;
		segByTarget.clear();

		if (path.empty())
			return false;
		const auto filePath = std::filesystem::path(std::u8string(path.begin(), path.end()));
		if (!std::filesystem::exists(filePath))
			return false;

		std::error_code sizeError;
		const auto fileSize = std::filesystem::file_size(filePath, sizeError);
		source = !sizeError && fileSize >= 16ull * 1024 * 1024 ? DataSourceFactory::openMmapped(path)
		                                                       : DataSourceFactory::open(path);
		if (source && !source->isOpen())
			source = DataSourceFactory::openLoaded(path);
		if (!source || !source->isOpen()) {
			source.reset();
			return false;
		}

		const auto bytes = source->bytes();
		if (bytes.size() < DAF_RECORD_BYTES) {
			LOG_WARN << "SpkReader: " << path << " — too small to be a DAF file (" << bytes.size() << " bytes)";
			source.reset();
			return false;
		}

		if (std::memcmp(bytes.data(), "DAF/SPK ", 8) != 0) {
			LOG_WARN << "SpkReader: " << path << " — not a DAF/SPK file (bad magic)";
			source.reset();
			return false;
		}

		byteSwap = (std::memcmp(bytes.data() + 88, "BIG-IEEE", 8) == 0);

		std::int32_t nd = 0, ni = 0, fward = 0;
		if (!readInt32At(bytes, 8, nd, byteSwap) || !readInt32At(bytes, 12, ni, byteSwap) ||
		    !readInt32At(bytes, 76, fward, byteSwap)) {
			LOG_WARN << "SpkReader: " << path << " — short read on file record";
			source.reset();
			return false;
		}

		if (nd != 2 || ni != 6) {
			LOG_WARN << "SpkReader: " << path << " — unexpected ND/NI (expected 2/6, got " << nd << "/" << ni << ")";
			source.reset();
			return false;
		}

		if (fward < 1) {
			LOG_WARN << "SpkReader: " << path << " — invalid FWARD=" << fward;
			source.reset();
			return false;
		}

		std::int32_t currRec = fward;
		while (currRec > 0) {
			const std::size_t recByteOffset = static_cast<std::size_t>(currRec - 1) * DAF_RECORD_BYTES;

			double header[3];
			if (!readDoublesAt(bytes, recByteOffset, 3, header, byteSwap))
				break;

			const std::int32_t nextRec = static_cast<std::int32_t>(header[0]);
			const std::int32_t nsum = static_cast<std::int32_t>(header[2]);

			std::size_t summaryOffset = recByteOffset + 3 * sizeof(double);

			bool chainTruncated = false;

			for (int s = 0; s < nsum; ++s) {
				double payload[SUMMARY_DOUBLES];
				if (!readDoublesAt(bytes, summaryOffset, SUMMARY_DOUBLES, payload, byteSwap)) {
					LOG_WARN << "SpkReader: unexpected EOF in summary record " << currRec;
					chainTruncated = true;
					break;
				}
				summaryOffset += SUMMARY_BYTES;

				Segment seg;
				seg.etStart = payload[0];
				seg.etEnd = payload[1];
				seg.target = doubleLowInt(payload[2], byteSwap);
				seg.center = doubleHighInt(payload[2], byteSwap);
				seg.frame = doubleLowInt(payload[3], byteSwap);
				seg.dataType = doubleHighInt(payload[3], byteSwap);
				seg.startAddr = doubleLowInt(payload[4], byteSwap);
				seg.endAddr = doubleHighInt(payload[4], byteSwap);

				static int s_diagCount = 0;
				if (s_diagCount < 20) {
					LOG_DEBUG << "SpkReader segment: target=" << seg.target << " center=" << seg.center
					          << " frame=" << seg.frame << " type=" << seg.dataType << " addr=[" << seg.startAddr
					          << ".." << seg.endAddr << "]";
					++s_diagCount;
				}

				if (seg.dataType != SPK_TYPE_CHEBYSHEV_POSITION)
					continue;
				if (seg.center != NAIF_SUN_ID && seg.center != NAIF_SSB_ID)
					continue;

				segByTarget[seg.target].push_back(seg);
			}

			if (chainTruncated)
				break;

			currRec = nextRec;
		}

		if (segByTarget.empty()) {
			LOG_WARN << "SpkReader: " << path << " — no Type-2 Sun/SSB-centered segments found";
			source.reset();
			return false;
		}

		std::size_t totalSegs = 0;
		int nSun = 0, nSsb = 0;
		for (auto& [tid, vec] : segByTarget) {
			std::sort(vec.begin(), vec.end(), [](const Segment& a, const Segment& b) { return a.etStart < b.etStart; });
			totalSegs += vec.size();
			for (const auto& seg : vec) {
				if (seg.center == NAIF_SUN_ID)
					++nSun;
				else if (seg.center == NAIF_SSB_ID)
					++nSsb;
			}
		}

		available = true;
		LOG_INFO << "SpkReader: loaded " << segByTarget.size() << " asteroids (" << totalSegs << " sub-segments) from "
		         << path << " [sun=" << nSun << ", ssb=" << nSsb << "]" << (byteSwap ? " [byte-swap]" : "");
		return true;
	}

	bool ensureTail(Segment& seg) const
	{
		if (seg.tailRead)
			return true;

		const auto bytes = view();
		if (bytes.empty()) {
			LOG_WARN << "SpkReader::ensureTail: source closed for " << sourcePath << " (target=" << seg.target << ")";
			return false;
		}

		const std::size_t tailOffset = static_cast<std::size_t>(seg.endAddr - 4) * sizeof(double);

		double tail[4];
		if (!readDoublesAt(bytes, tailOffset, 4, tail, byteSwap)) {
			LOG_WARN << "SpkReader::ensureTail: read of 4 doubles at offset " << tailOffset
			         << " failed (target=" << seg.target << " endAddr=" << seg.endAddr << " view.size=" << bytes.size()
			         << ")";
			return false;
		}

		seg.initEt = tail[0];
		seg.intLen = tail[1];
		seg.rsize = static_cast<int>(tail[2]);
		seg.nRecords = static_cast<int>(tail[3]);

		if (seg.intLen <= 0.0 || seg.rsize < 5 || seg.nRecords < 1) {
			LOG_WARN << "SpkReader: target " << seg.target << " has invalid Type-2 tail (INTLEN=" << seg.intLen
			         << ", RSIZE=" << seg.rsize << ", N=" << seg.nRecords << ")";
			return false;
		}
		if ((seg.rsize - 2) % 3 != 0) {
			LOG_WARN << "SpkReader: target " << seg.target << " Type-2 RSIZE=" << seg.rsize << " not of form 2+3*NC";
			return false;
		}

		seg.nc = (seg.rsize - 2) / 3;
		seg.tailRead = true;
		return true;
	}

	bool evaluate(const Segment& segIn, orbital_core::detail::SplitTime et, Vec3D& R_km, Vec3D* V_km_s) const
	{
		Segment& seg = const_cast<Segment&>(segIn);
		if (!ensureTail(seg))
			return false;

		if (!et.atOrAfter(seg.etStart) || !et.atOrBefore(seg.etEnd))
			return false;

		int idx = static_cast<int>(std::floor(et.since(seg.initEt) / seg.intLen));
		if (idx < 0)
			idx = 0;
		if (idx >= seg.nRecords)
			idx = seg.nRecords - 1;
		if (idx > 0 && !et.atOrAfter(seg.initEt + idx * seg.intLen))
			--idx;

		const std::size_t recOffset = static_cast<std::size_t>(seg.startAddr - 1 + idx * seg.rsize) * sizeof(double);

		const auto bytes = view();
		if (bytes.empty()) {
			LOG_WARN << "SpkReader::evaluate: source closed (target=" << seg.target << ")";
			return false;
		}

		constexpr int kMaxStackRSize = 64;
		std::array<double, kMaxStackRSize> stackBuf;
		std::vector<double> heapBuf;
		const std::size_t recSize = static_cast<std::size_t>(seg.rsize);

		double* recData = nullptr;
		if (seg.rsize <= kMaxStackRSize)
			recData = stackBuf.data();

		else {
			heapBuf.resize(recSize);
			recData = heapBuf.data();
		}

		if (!readDoublesAt(bytes, recOffset, recSize, recData, byteSwap)) {
			LOG_WARN << "SpkReader::evaluate: read of " << seg.rsize << " doubles at offset " << recOffset
			         << " failed (target=" << seg.target << " idx=" << idx << " view.size=" << bytes.size() << ")";
			return false;
		}

		const double mid = recData[0];
		const double radius = recData[1];
		if (radius <= 0.0)
			return false;

		const int nc = seg.nc;
		const double* CX = recData + 2;
		const double* CY = CX + nc;
		const double* CZ = CY + nc;

		const double tau = et.since(mid) / radius;

		auto evalAxis = [&](const double* C, double& pos, double& velTau) {
			double b1 = 0.0, b2 = 0.0;
			double bp1 = 0.0, bp2 = 0.0;
			const double twoTau = 2.0 * tau;
			for (int k = nc - 1; k >= 1; --k) {
				const double bk_new = twoTau * b1 - b2 + C[k];
				const double bpk_new = twoTau * bp1 - bp2 + 2.0 * b1;
				b2 = b1;
				b1 = bk_new;
				bp2 = bp1;
				bp1 = bpk_new;
			}

			pos = tau * b1 - b2 + C[0];
			velTau = tau * bp1 - bp2 + b1;
		};

		double px, py, pz, vTauX, vTauY, vTauZ;
		evalAxis(CX, px, vTauX);
		evalAxis(CY, py, vTauY);
		evalAxis(CZ, pz, vTauZ);

		R_km = Vec3D(px, py, pz);

		if (V_km_s) {
			*V_km_s = Vec3D(vTauX / radius, vTauY / radius, vTauZ / radius);
		}
		return true;
	}
};

SpkReader::SpkReader() : d(std::make_unique<Impl>())
{
	static const char* candidates[] = {"sb441-n16.bsp", "sb441-n343.bsp"};

	for (const auto& directory : {orbital_core::ephemerisDataDirectory(), getExecutableDir()}) {
		if (directory.empty())
			continue;
		for (const char* fn : candidates) {
			const auto utf8 = (directory / fn).u8string();
			if (d->open(std::string(utf8.begin(), utf8.end())))
				return;
		}
	}
}

SpkReader::SpkReader(const std::string& bspPath) : d(std::make_unique<Impl>()) { d->open(bspPath); }

SpkReader::~SpkReader() = default;

bool SpkReader::isAvailable() const noexcept { return d->available; }

std::string SpkReader::sourcePath() const { return d->available ? d->sourcePath : std::string{}; }

std::vector<int> SpkReader::targetIds() const
{
	std::vector<int> ids;
	ids.reserve(d->segByTarget.size());

	for (const auto& [tid, _] : d->segByTarget)
		ids.push_back(tid);

	return ids;
}

static Vec3D sunFromSsbEclAu(double) { return {0.0, 0.0, 0.0}; }

bool SpkReader::queryState(int naifId, double jdTdb, Vec3D& R, Vec3D& V) const
{
	return queryStateMjd(naifId, jdTdb - 2400000.5, R, V);
}

bool SpkReader::queryStateMjd(int naifId, double mjdTdb, Vec3D& R, Vec3D& V) const
{
	return queryStateOffsetMjd(naifId, mjdTdb, 0.0, R, V);
}

bool SpkReader::queryStateOffsetMjd(int naifId, double epochMjd, double offsetDays, Vec3D& R, Vec3D& V) const
{
	if (!d->available)
		return false;

	const auto et = orbital_core::detail::SplitTime::sum(epochMjd, offsetDays).shifted(51544.5).scaled(86400.0);
	const auto* seg = d->pickSegment(naifId, et);
	if (!seg)
		return false;

	Vec3D R_km, V_km_s;
	if (!d->evaluate(*seg, et, R_km, &V_km_s))
		return false;

	R = rotateEqToEcl(R_km.x / AU_KM, R_km.y / AU_KM, R_km.z / AU_KM);

	const double kmsToAuDay = 86400.0 / AU_KM;
	V = rotateEqToEcl(V_km_s.x * kmsToAuDay, V_km_s.y * kmsToAuDay, V_km_s.z * kmsToAuDay);

	if (seg->center == NAIF_SSB_ID)
		R = R - sunFromSsbEclAu(epochMjd + offsetDays + 2400000.5);

	return true;
}

Vec3D SpkReader::getPosition(int naifId, double jdTdb) const { return getPositionMjd(naifId, jdTdb - 2400000.5); }

Vec3D SpkReader::getPositionMjd(int naifId, double mjdTdb) const { return getPositionOffsetMjd(naifId, mjdTdb, 0.0); }

Vec3D SpkReader::getPositionOffsetMjd(int naifId, double epochMjd, double offsetDays) const
{
	if (!d->available)
		return {0.0, 0.0, 0.0};

	const auto et = orbital_core::detail::SplitTime::sum(epochMjd, offsetDays).shifted(51544.5).scaled(86400.0);
	const auto* seg = d->pickSegment(naifId, et);
	if (!seg)
		return {0.0, 0.0, 0.0};

	Vec3D R_km;
	if (!d->evaluate(*seg, et, R_km, nullptr))
		return {0.0, 0.0, 0.0};

	Vec3D R = rotateEqToEcl(R_km.x / AU_KM, R_km.y / AU_KM, R_km.z / AU_KM);

	if (seg->center == NAIF_SSB_ID)
		R = R - sunFromSsbEclAu(epochMjd + offsetDays + 2400000.5);

	return R;
}

static SpkReader* s_defaultAsteroidSource = nullptr;

SpkReader* getDefaultAsteroidSource() { return s_defaultAsteroidSource; }

void setDefaultAsteroidSource(SpkReader* source) { s_defaultAsteroidSource = source; }
