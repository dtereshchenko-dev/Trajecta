#include "Sb441Download.h"

#include <QDir>
#include <QFileInfo>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QSaveFile>
#include <QStandardPaths>
#include <algorithm>
#include <utility>

namespace {
	constexpr qint64 kSb441Bytes = 645727232;
	const auto kSb441Sha256 = QByteArray::fromHex("919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90");
	const QUrl kSb441Url(QStringLiteral("https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp"));

	QString installedDirectory()
	{
		const auto root = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
		return root.isEmpty() ? QString() : QDir(root).filePath("ephemerides");
	}
}  // namespace

QString Sb441Download::installedPath()
{
	const auto directory = installedDirectory();
	return directory.isEmpty() ? QString() : QDir(directory).filePath("sb441-n16.bsp");
}

Sb441Download::Sb441Download(QObject* parent)
    : Sb441Download(kSb441Url, installedDirectory(), kSb441Bytes, kSb441Sha256, parent)
{
}

Sb441Download::Sb441Download(QUrl source, QString directory, qint64 expectedSize, QByteArray expectedSha256,
                             QObject* parent)
    : QObject(parent),
      m_source(std::move(source)),
      m_directory(std::move(directory)),
      m_expectedSize(expectedSize),
      m_expectedSha256(std::move(expectedSha256))
{
}

Sb441Download::~Sb441Download()
{
	if (m_reply) {
		m_reply->disconnect(this);
		m_reply->abort();
	}
	if (m_file)
		m_file->cancelWriting();
}

double Sb441Download::progress() const
{
	return m_expectedSize > 0 ? std::min(1.0, static_cast<double>(m_received) / m_expectedSize) : 0.0;
}

void Sb441Download::start()
{
	if (running())
		return;
	m_cancelled = false;
	m_failure.clear();
	m_received = 0;
	m_hash.reset();
	if (m_directory.isEmpty() || !QDir().mkpath(m_directory)) {
		m_message = tr("Cannot create the data directory.");
		emit changed();
		return;
	}
	m_file = std::make_unique<QSaveFile>(QDir(m_directory).filePath("sb441-n16.bsp"));
	m_file->setDirectWriteFallback(false);
	if (!m_file->open(QIODevice::WriteOnly)) {
		m_message = tr("Cannot create the temporary download: %1").arg(m_file->errorString());
		m_file.reset();
		emit changed();
		return;
	}
	QNetworkRequest request(m_source);
	request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::NoLessSafeRedirectPolicy);
	request.setTransferTimeout(60000);
	m_reply = m_network.get(request);
	m_reply->setReadBufferSize(1024 * 1024);
	connect(m_reply, &QIODevice::readyRead, this, &Sb441Download::consume);
	connect(m_reply, &QNetworkReply::finished, this, &Sb441Download::finish);
	m_message = tr("Downloading SB441-N16...");
	emit changed();
}

void Sb441Download::cancel()
{
	if (!m_reply)
		return;
	m_cancelled = true;
	m_reply->abort();
}

void Sb441Download::consume()
{
	if (!m_reply || !m_file || m_cancelled || !m_failure.isEmpty())
		return;
	while (m_reply->bytesAvailable() > 0) {
		const auto chunk = m_reply->read(256 * 1024);
		if (chunk.isEmpty())
			break;
		if (m_received > m_expectedSize - chunk.size()) {
			m_failure = tr("The download is larger than the expected JPL file.");
			if (!m_finishing)
				m_reply->abort();
			return;
		}
		if (m_file->write(chunk) != chunk.size()) {
			m_failure = tr("Cannot write the download: %1").arg(m_file->errorString());
			if (!m_finishing)
				m_reply->abort();
			return;
		}
		m_hash.addData(chunk);
		m_received += chunk.size();
	}
	m_message = tr("Downloading: %1 / %2 MiB")
	                .arg(m_received / 1048576.0, 0, 'f', 1)
	                .arg(m_expectedSize / 1048576.0, 0, 'f', 1);
	emit changed();
}

void Sb441Download::finish()
{
	if (!m_reply)
		return;
	m_finishing = true;
	consume();
	const auto status = m_reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
	const auto finalUrl = m_reply->url();
	if (m_cancelled)
		m_message = tr("Download cancelled.");
	else if (!m_failure.isEmpty())
		m_message = m_failure;
	else if (m_reply->error() != QNetworkReply::NoError)
		m_message = tr("Download failed: %1").arg(m_reply->errorString());
	else if (status != 200 || finalUrl.scheme() != m_source.scheme() || finalUrl.host() != m_source.host())
		m_message = tr("Unexpected response or download destination.");
	else if (m_received != m_expectedSize || m_hash.result() != m_expectedSha256)
		m_message = tr("SB441 verification failed. The existing file was not changed.");
	else if (!m_file->commit())
		m_message = tr("Cannot install SB441: %1").arg(m_file->errorString());
	else {
		m_message = tr("SB441 downloaded and verified. Restart the application to use it.");
		m_file.reset();
		m_reply->deleteLater();
		m_reply = nullptr;
		m_finishing = false;
		emit changed();
		emit installed();
		return;
	}
	if (m_file)
		m_file->cancelWriting();
	m_file.reset();
	m_reply->deleteLater();
	m_reply = nullptr;
	m_finishing = false;
	emit changed();
}
