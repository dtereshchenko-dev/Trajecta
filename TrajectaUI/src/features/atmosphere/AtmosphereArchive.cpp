#include "AtmosphereArchive.h"

#include <QCryptographicHash>
#include <QDateTime>
#include <QDirIterator>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QSaveFile>
#include <QTimeZone>
#include <orbital_core/OrbitSolution.h>
#include <orbital_core/TimeScale.h>
#include <orbital_core/hazard/BodyFixedFrame.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace atmosphere {
	namespace {
		void require(bool condition, const char* message)
		{
			if (!condition)
				throw std::invalid_argument(message);
		}

		QByteArray read(const QString& path)
		{
			QFile file(path);
			if (!file.open(QIODevice::ReadOnly))
				throw std::runtime_error(("Cannot read " + path).toStdString());
			const auto bytes = file.readAll();
			if (file.error() != QFile::NoError)
				throw std::runtime_error(("Cannot read " + path).toStdString());
			return bytes;
		}

		QByteArray hash(const QByteArray& data)
		{
			return QCryptographicHash::hash(data, QCryptographicHash::Sha256).toHex();
		}

		double number(const QJsonObject& object, const char* key)
		{
			const auto value = object.value(QLatin1String(key));
			require(value.isDouble() && std::isfinite(value.toDouble()), "Missing or non-finite archive number");
			return value.toDouble();
		}

		int count(const QJsonObject& object, const char* key)
		{
			const auto value = number(object, key);
			require(value >= 0 && value <= 100000000 && std::floor(value) == value, "Invalid archive count");
			return static_cast<int>(value);
		}

		Vec3D vector(const QJsonValue& value)
		{
			require(value.isArray() && value.toArray().size() == 3, "Invalid entry vector");
			const auto a = value.toArray();
			for (const auto& v : a)
				require(v.isDouble() && std::isfinite(v.toDouble()), "Non-finite entry vector");
			return {a[0].toDouble(), a[1].toDouble(), a[2].toDouble()};
		}

		strewn::Vector3 array(const Vec3D& v) { return {v.x, v.y, v.z}; }

		QJsonArray json(const strewn::Vector3& v) { return {v[0], v[1], v[2]}; }
	}  // namespace

	double parseUtc(const QString& text)
	{
		const QDateTime date(QDate::fromString(text.left(10), QStringLiteral("yyyy-MM-dd")),
		                     QTime::fromString(text.mid(11), QStringLiteral("HH:mm:ss")), QTimeZone::UTC);
		require(date.isValid() && date.toString("yyyy-MM-dd HH:mm:ss") == text, "Enter YYYY-MM-DD HH:MM:SS UTC");
		const double utc = 40587.0 + date.toMSecsSinceEpoch() / 86400000.0;
		require(std::isfinite(orbital_core::time::toMjdTdb(utc, orbital_core::time::Scale::UTC)),
		        "UTC epoch is unavailable");
		return utc;
	}

	QString formatTdb(double tdb)
	{
		const double utc = orbital_core::time::fromMjdTdb(tdb, orbital_core::time::Scale::UTC);
		return QDateTime::fromMSecsSinceEpoch(qRound64((utc - 40587.0) * 86400000.0), QTimeZone::UTC)
		    .toString(QStringLiteral("yyyy-MM-dd HH:mm:ss.zzz"));
	}

	void validateEntry(const strewn::EntryState& e)
	{
		for (const auto& a : {e.positionEcefM, e.velocityEcefMps})
			for (double v : a)
				require(std::isfinite(v), "Non-finite entry state");
		require(std::isfinite(e.epochMjdUtc) && std::isfinite(e.weight) && e.weight > 0,
		        "Invalid entry time or weight");
		const auto geo = strewn::geodetic(e.positionEcefM);
		require(std::isfinite(geo.heightWgs84M) && std::abs(geo.heightWgs84M - 100000) <= 1.0,
		        "Entry is not at WGS84 height 100 km");
		const auto up = strewn::enuToEcef({0, 0, 1}, geo.latitudeRad, geo.longitudeRad);
		double vertical = 0;
		for (int i = 0; i < 3; ++i)
			vertical += up[i] * e.velocityEcefMps[i];
		require(vertical < 0, "Entry is not inbound");
	}

	EntryArchive readArchive(const QString& path, std::stop_token stop)
	{
		namespace hz = orbital_core::hazard;
		EntryArchive result;
		result.path = QFileInfo(path).absoluteFilePath();
		const auto bytes = read(path);
		result.hash = hash(bytes);
		QJsonParseError parseError;
		const auto document = QJsonDocument::fromJson(bytes, &parseError);
		require(parseError.error == QJsonParseError::NoError && document.isObject(), "Invalid MC JSON archive");
		const auto root = document.object();
		const bool tdbArchive = root["schema"] == "orbital-interval-mc-4";
		require(tdbArchive || root["schema"] == "orbital-interval-mc-3",
		        "Expected an interval MC version 3 or 4 archive");
		require(root["time_scale"] == (tdbArchive ? "MJD TDB" : "MJD TT") &&
		            root["frame"] == "heliocentric ecliptic J2000" && root["state_units"] == "AU, AU/day",
		        "Unsupported MC time scale, frame or units");
		const auto epochTdb = [tdbArchive](double value) {
			return orbital_core::time::toMjdTdb(
			    value, tdbArchive ? orbital_core::time::Scale::TDB : orbital_core::time::Scale::TT);
		};
		require(count(root, "target") == static_cast<int>(hz::Planet::Earth) && number(root, "altitude_km") == 100,
		        "Atmospheric MC requires Earth entries at 100 km");
		require(root["body_fixed_frame"] == "EARTH_OPTICAL_UT1" &&
		            root["body_fixed_model"].toString().toStdString() == hz::earthBodyFixedFrameModel(),
		        "Archive Earth orientation model differs from the active model");
		const auto requested = count(root, "requested");
		require(root["cancelled"].isBool() && !root["cancelled"].toBool() && requested > 0 &&
		            count(root, "completed") == requested && count(root, "failed_fit") == 0 &&
		            count(root, "failed_propagation") == 0 && count(root, "cancelled_clones") == 0 &&
		            count(root, "unstarted") == 0,
		        "Atmospheric MC requires a complete successful orbital ensemble");
		require(root["clones"].isArray() && root["clones"].toArray().size() == requested && root["events"].isArray(),
		        "Incomplete MC archive");
		require(root["method"] == "covariance" || root["method"] == "observations", "Unknown MC sampling method");
		std::vector<int> expected;
		const double start = epochTdb(number(root, "start")), end = epochTdb(number(root, "end"));
		for (const auto& value : root["events"].toArray()) {
			const auto event = value.toObject();
			EntryEvent e{epochTdb(number(event, "start")), epochTdb(number(event, "end"))};
			require(e.startMjdTdb <= e.endMjdTdb && e.startMjdTdb >= start && e.endMjdTdb <= end &&
			            (result.events.empty() || e.startMjdTdb > result.events.back().endMjdTdb),
			        "Invalid or overlapping MC events");
			expected.push_back(count(event, "entries"));
			require(expected.back() > 0, "Empty MC event");
			result.events.push_back(std::move(e));
		}
		std::set<int> ids;
		int entries = 0;
		for (const auto& value : root["clones"].toArray()) {
			if (stop.stop_requested())
				throw std::runtime_error("Archive loading cancelled");
			const auto row = value.toObject();
			const int id = count(row, "index");
			require(id < requested && ids.insert(id).second && count(row, "status") == 1,
			        "Duplicate or incomplete MC clone");
			if (row["entry"].isNull())
				continue;
			require(row["entry"].isObject() && row["entry_before_interval"].isBool() &&
			            !row["entry_before_interval"].toBool(),
			        "Invalid crossing record");
			const auto entry = row["entry"].toObject();
			const double epoch = epochTdb(number(entry, "epoch"));
			auto event = std::find_if(result.events.begin(), result.events.end(), [epoch](const auto& e) {
				return epoch >= e.startMjdTdb && epoch <= e.endMjdTdb;
			});
			require(event != result.events.end(), "Crossing is outside the saved event intervals");
			const auto r = vector(entry["relative_position_km"]), v = vector(entry["relative_velocity_km_s"]);
			const auto position = hz::icrfToBodyFixedMjd(hz::Planet::Earth, r, epoch);
			const double before = epoch - 1.0 / 86400, after = epoch + 1.0 / 86400;
			const auto rate = (hz::icrfToBodyFixedMjd(hz::Planet::Earth, r, after) -
			                   hz::icrfToBodyFixedMjd(hz::Planet::Earth, r, before)) /
			                  ((after - before) * 86400);
			const auto velocity = hz::icrfToBodyFixedMjd(hz::Planet::Earth, v, epoch) + rate;
			strewn::EntryState e{array(position * 1000), array(velocity * 1000),
			                     orbital_core::time::fromMjdTdb(epoch, orbital_core::time::Scale::UTC), 1.0,
			                     static_cast<std::uint64_t>(id)};
			validateEntry(e);
			const auto geo = strewn::geodetic(e.positionEcefM);
			constexpr double degrees = 180 / 3.14159265358979323846;
			require(
			    std::abs(geo.latitudeRad * degrees - number(entry, "latitude_deg")) < 1e-6 &&
			        std::abs(std::remainder(geo.longitudeRad * degrees - number(entry, "longitude_deg"), 360)) < 1e-6,
			    "Entry coordinates differ from the saved Earth orientation; regenerate the orbital archive");
			event->entries.push_back(e);
			++entries;
		}
		require(entries > 0 && entries == count(root, "entry_clones"), "Inconsistent or empty crossing ensemble");
		for (std::size_t i = 0; i < result.events.size(); ++i)
			require(result.events[i].entries.size() == expected[i], "MC event crossing count mismatch");
		const auto work = QFileInfo(path).absoluteDir().filePath("work");
		QDirIterator oefs(work, {"orbit_equinoctial_oef.txt"}, QDir::Files, QDirIterator::Subdirectories);
		if (oefs.hasNext()) {
			result.photometryPath = oefs.next();
			require(!oefs.hasNext(), "Ambiguous source solution inside MC archive");
			result.photometryHash = hash(read(result.photometryPath));
			result.magnitude =
			    orbital_core::loadOrbitSolution(std::filesystem::path(result.photometryPath.toStdWString()))
			        .absoluteMagnitude;
		}
		result.provenance = {
		    {"schema", "atmospheric-interval-handoff-1"},
		    {"source_path", result.path},
		    {"source_sha256", QString::fromLatin1(result.hash)},
		    {"orbital_method", root["method"]},
		    {"orbital_requested", requested},
		    {"orbital_completed", requested},
		    {"orbital_entries", entries},
		    {"source_time_scale", root["time_scale"]},
		    {"source_entry_frame", "geocentric ecliptic J2000"},
		    {"source_entry_units", "km, km/s"},
		    {"entry_frame", "ECEF"},
		    {"entry_units", "m, m/s"},
		    {"earth_orientation", root["body_fixed_model"]},
		    {"velocity_transform", "R*v + dR/dt*r; symmetric 1-second frame derivative"},
		    {"conditioning", "selected event entry clones; not an unconditional Earth impact probability"},
		    {"photometry_path", result.photometryPath},
		    {"photometry_sha256", QString::fromLatin1(result.photometryHash)}};
		return result;
	}

	void verifyArchive(const EntryArchive& archive)
	{
		require(hash(read(archive.path)) == archive.hash, "MC archive changed after selection; select it again");
		if (!archive.photometryPath.isEmpty())
			require(hash(read(archive.photometryPath)) == archive.photometryHash,
			        "Source photometry changed; select the archive again");
		require(archive.provenance["earth_orientation"].toString().toStdString() ==
		            orbital_core::hazard::earthBodyFixedFrameModel(),
		        "Earth orientation changed after archive selection");
	}

	void saveJson(const QString& path, const QJsonObject& value)
	{
		QSaveFile file(path);
		if (!file.open(QIODevice::WriteOnly))
			throw std::runtime_error(("Cannot write " + path).toStdString());
		const auto data = QJsonDocument(value).toJson();
		if (file.write(data) != data.size() || !file.commit())
			throw std::runtime_error(("Cannot save " + path).toStdString());
	}

	void saveEntries(const QString& path, const std::vector<strewn::EntryState>& entries)
	{
		QJsonArray rows;
		for (const auto& e : entries)
			rows.push_back(QJsonObject{{"id", qint64(e.id)},
			                           {"epoch_mjd_utc", e.epochMjdUtc},
			                           {"position_ecef_m", json(e.positionEcefM)},
			                           {"velocity_ecef_mps", json(e.velocityEcefMps)},
			                           {"weight", e.weight}});
		saveJson(path, {{"schema_version", 1}, {"frame", "ECEF"}, {"altitude_km", 100}, {"entries", rows}});
	}
}  // namespace atmosphere
