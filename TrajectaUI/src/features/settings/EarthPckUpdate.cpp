#include "EarthPckUpdate.h"

#include "orbital_core/RuntimeData.h"

#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSettings>
#include <QStandardPaths>
#include <QTemporaryFile>
#include <algorithm>

namespace {
	constexpr qint64 maxPckBytes = 128 * 1024 * 1024;
	const QUrl naifDirectory("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/");
	const std::array<std::pair<const char*, const char*>, 3> pckFiles{{{"earth_pck_historic", "earth_historic.bpc"},
	                                                                   {"earth_pck_latest", "earth_latest.bpc"},
	                                                                   {"earth_pck_predict", "earth_predict.bpc"}}};

	QString pckDirectory() { return QString::fromStdWString(orbital_core::earthPckDataDirectory().wstring()); }

	QString metadataPath()
	{
		return QDir(QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation))
		    .filePath("earth-pck-updates.ini");
	}

	QString fileHash(const QString& path)
	{
		QFile file(path);
		if (!file.open(QIODevice::ReadOnly))
			return {};
		QCryptographicHash hash(QCryptographicHash::Sha256);
		return hash.addData(&file) ? QString::fromLatin1(hash.result().toHex()) : QString();
	}

	QString remoteName(const QString& listing, const QString& pattern)
	{
		QRegularExpression re(pattern);
		auto matches = re.globalMatch(listing);
		QString result;
		while (matches.hasNext()) {
			const auto candidate = matches.next().captured(1);
			if (candidate > result)
				result = candidate;
		}
		return result;
	}

	bool responseOk(QNetworkReply* reply, const QUrl& source)
	{
		return reply->error() == QNetworkReply::NoError &&
		       reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt() == 200 &&
		       reply->url().scheme() == source.scheme() && reply->url().host() == source.host();
	}
}  // namespace

EarthPckUpdate::EarthPckUpdate(QObject* parent) : EarthPckUpdate(naifDirectory, pckDirectory(), metadataPath(), parent)
{
}

void EarthPckUpdate::applyPending() { applyPending(pckDirectory(), metadataPath()); }

void EarthPckUpdate::applyPending(const QString& directory, const QString& metadataPath)
{
	QSettings metadata(metadataPath, QSettings::IniFormat);
	for (const auto& [id, file] : pckFiles) {
		const QString key = QString::fromLatin1(id);
		const auto name = metadata.value(key + "/pendingName").toString();
		if (name.isEmpty())
			continue;
		const auto destination = QDir(directory).filePath(QString::fromLatin1(file));
		const auto pending = destination + ".pending";
		const auto expected = metadata.value(key + "/pendingHash").toString();
		if (expected.isEmpty() || fileHash(pending) != expected ||
		    fileHash(destination) != metadata.value(key + "/pendingOriginalHash").toString()) {
			metadata.setValue(key + "/pendingError", "Pending update or local PCK changed; check again.");
			continue;
		}
		QFile input(pending);
		QSaveFile output(destination);
		output.setDirectWriteFallback(false);
		if (!input.open(QIODevice::ReadOnly) || !output.open(QIODevice::WriteOnly)) {
			metadata.setValue(key + "/pendingError", "Cannot apply the pending PCK update.");
			continue;
		}
		bool ok = true;
		QCryptographicHash hash(QCryptographicHash::Sha256);
		while (!input.atEnd()) {
			const auto block = input.read(256 * 1024);
			if (block.isEmpty() || output.write(block) != block.size()) {
				ok = false;
				break;
			}
			hash.addData(block);
		}
		if (!ok || input.error() != QFileDevice::NoError || QString::fromLatin1(hash.result().toHex()) != expected ||
		    fileHash(destination) != metadata.value(key + "/pendingOriginalHash").toString() || !output.commit()) {
			metadata.setValue(key + "/pendingError", "Cannot apply the pending PCK update.");
			continue;
		}
		input.close();
		QFile::remove(pending);
		metadata.setValue(key + "/name", name);
		metadata.setValue(key + "/sha256", expected);
		metadata.remove(key + "/pendingName");
		metadata.remove(key + "/pendingHash");
		metadata.remove(key + "/pendingOriginalHash");
		metadata.remove(key + "/pendingError");
	}
	metadata.sync();
}

EarthPckUpdate::EarthPckUpdate(QUrl directoryUrl, QString directory, QString metadataPath, QObject* parent)
    : QObject(parent),
      m_directoryUrl(std::move(directoryUrl)),
      m_directory(std::move(directory)),
      m_metadataPath(std::move(metadataPath)),
      m_sources{{{"earth_pck_historic", "earth_historic.bpc"},
                 {"earth_pck_latest", "earth_latest.bpc"},
                 {"earth_pck_predict", "earth_predict.bpc"}}}
{
}

EarthPckUpdate::Source* EarthPckUpdate::source(const QString& id)
{
	for (auto& item : m_sources)
		if (item.id == id)
			return &item;
	return nullptr;
}

QVariantMap EarthPckUpdate::states() const
{
	QVariantMap result;
	for (const auto& item : m_sources)
		result.insert(item.id, QVariantMap{{"available", item.available},
		                                   {"checking", item.checking},
		                                   {"downloading", item.downloading},
		                                   {"message", item.message}});
	return result;
}

int EarthPckUpdate::availableCount() const
{
	return static_cast<int>(
	    std::count_if(m_sources.begin(), m_sources.end(), [](const auto& s) { return s.available; }));
}

QNetworkReply* EarthPckUpdate::get(const QUrl& url)
{
	QNetworkRequest request(url);
	request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::NoLessSafeRedirectPolicy);
	request.setTransferTimeout(60000);
	return m_network.get(request);
}

void EarthPckUpdate::check()
{
	if (m_check)
		return;
	for (auto& item : m_sources) {
		if (!item.downloading) {
			item.checking = true;
			item.available = false;
			item.message.clear();
		}
	}
	m_check = get(m_directoryUrl);
	connect(m_check, &QNetworkReply::finished, this, [this, reply = m_check.data()] { finishCheck(reply); });
	emit changed();
}

void EarthPckUpdate::finishCheck(QNetworkReply* reply)
{
	const bool ok = responseOk(reply, m_directoryUrl) && reply->bytesAvailable() <= 2 * 1024 * 1024;
	const QString listing = ok ? QString::fromUtf8(reply->readAll()) : QString();
	const std::array<QString, 3> names{remoteName(listing, R"((earth_\d{6}_\d{6}\.bpc))"),
	                                   remoteName(listing, R"((earth_000101_\d{6}_\d{6}\.bpc))"),
	                                   remoteName(listing, R"((earth_\d{4}_\d{6}_\d{4}_predict\.bpc))")};
	QSettings metadata(m_metadataPath, QSettings::IniFormat);
	const std::array<QString, 3> initialNames{"earth_620120_260806.bpc", "earth_000101_270208_261112.bpc",
	                                          "earth_2026_260806_2126_predict.bpc"};
	const std::array<QString, 3> initialHashes{"652b02c0b3a5a283d393497f7f96e1fe4aa47ee6dcc035a089eb9e6b5df86643",
	                                           "eb72aef57bcedb247e282510a881a608fe74e4d9e7869aaa130474f5efe8dc54",
	                                           "9f565ab5aaccd33afe54f3a1b20cb8b136636c31a0bd4c0986fac722b931e280"};
	for (std::size_t i = 0; i < m_sources.size(); ++i) {
		auto& item = m_sources[i];
		if (item.downloading)
			continue;
		item.checking = false;
		if (!metadata.value(item.id + "/pendingName").toString().isEmpty() &&
		    metadata.value(item.id + "/pendingError").toString().isEmpty()) {
			item.message = tr("Update ready. Restart the application to apply it.");
			continue;
		}
		if (!ok || names[i].isEmpty()) {
			item.message = tr("Update check unavailable.");
			continue;
		}
		item.remote = names[i];
		const auto localPath = QDir(m_directory).filePath(item.file);
		item.originalHash = fileHash(localPath);
		const auto activePath =
		    m_directory == pckDirectory()
		        ? QString::fromStdWString(orbital_core::earthPckFilePath(item.file.toStdString()).wstring())
		        : localPath;
		const auto activeHash = fileHash(activePath);
		if ((QFileInfo::exists(localPath) && item.originalHash.isEmpty()) ||
		    (QFileInfo::exists(activePath) && activeHash.isEmpty())) {
			item.message = tr("Cannot read the local PCK file.");
			continue;
		}
		const bool savedMatch = metadata.value(item.id + "/name").toString() == item.remote &&
		                        metadata.value(item.id + "/sha256").toString() == activeHash;
		const bool initialMatch = item.remote == initialNames[i] && activeHash == initialHashes[i];
		item.available = activeHash.isEmpty() || (!savedMatch && !initialMatch) ||
		                 !metadata.value(item.id + "/pendingError").toString().isEmpty();
		if (item.available) {
			item.message = metadata.value(item.id + "/pendingError").toString();
			if (!item.message.isEmpty())
				item.message += QLatin1Char(' ');
			item.message += tr("Available: %1").arg(item.remote);
		}
	}
	reply->deleteLater();
	m_check = nullptr;
	emit changed();
}

void EarthPckUpdate::install(const QString& id)
{
	auto* item = source(id);
	if (!item || !item->available || item->checking || item->downloading)
		return;
	item->downloading = true;
	item->message = tr("Downloading %1...").arg(item->remote);
	auto* reply = get(m_directoryUrl.resolved(QUrl(item->remote)));
	connect(reply, &QNetworkReply::downloadProgress, this, [reply](qint64 received, qint64 total) {
		if (received > maxPckBytes || total > maxPckBytes)
			reply->abort();
	});
	connect(reply, &QNetworkReply::finished, this, [this, reply, item] { finishDownload(reply, item); });
	emit changed();
}

void EarthPckUpdate::finishDownload(QNetworkReply* reply, Source* item)
{
	const auto destination = QDir(m_directory).filePath(item->file);
	const auto data = reply->readAll();
	QString error;
	const auto contentLength = reply->header(QNetworkRequest::ContentLengthHeader);
	if (!responseOk(reply, m_directoryUrl) || (contentLength.isValid() && contentLength.toLongLong() != data.size()) ||
	    data.size() > maxPckBytes || data.size() < 1024 || !data.startsWith("DAF/PCK"))
		error = tr("Download failed or is not a binary PCK.");
	else if (fileHash(destination) != item->originalHash)
		error = tr("Local file changed; check for updates again.");
	else if (!QDir().mkpath(m_directory) || !QDir().mkpath(QFileInfo(m_metadataPath).absolutePath()))
		error = tr("Cannot create the Earth orientation or update metadata directory.");
	else {
		QTemporaryFile temp(QDir(m_directory).filePath(".pck-check-XXXXXX"));
		if (!temp.open() || temp.write(data) != data.size() || !temp.flush())
			error = tr("Cannot prepare the downloaded PCK.");
		else {
			const auto tempPath = temp.fileName();
			temp.close();
			if (!orbital_core::validateEarthPck(std::filesystem::path(tempPath.toStdWString())))
				error = tr("PCK validation failed; the existing file was not changed.");
			else {
				QSaveFile output(destination + ".pending");
				output.setDirectWriteFallback(false);
				if (!output.open(QIODevice::WriteOnly) || output.write(data) != data.size() || !output.commit())
					error = tr("Cannot prepare PCK update: %1").arg(output.errorString());
				else {
					QSettings metadata(m_metadataPath, QSettings::IniFormat);
					metadata.setValue(item->id + "/pendingName", item->remote);
					metadata.setValue(item->id + "/pendingOriginalHash", item->originalHash);
					metadata.setValue(
					    item->id + "/pendingHash",
					    QString::fromLatin1(QCryptographicHash::hash(data, QCryptographicHash::Sha256).toHex()));
					metadata.remove(item->id + "/pendingError");
					metadata.sync();
					if (metadata.status() != QSettings::NoError)
						error = tr("Cannot save PCK update metadata.");
					else {
						item->available = false;
						item->message = tr("Update ready. Restart the application to apply it.");
						emit installed();
					}
				}
			}
		}
	}
	if (!error.isEmpty())
		item->message = error;
	item->downloading = false;
	reply->deleteLater();
	emit changed();
}
