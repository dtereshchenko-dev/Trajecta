#include "ApplicationUpdate.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QVersionNumber>
#include <QXmlStreamReader>

ApplicationUpdate::ApplicationUpdate(QObject* parent) : QObject(parent)
{
	m_timeout.setSingleShot(true);
	connect(&m_timeout, &QTimer::timeout, this, [this] { m_process.kill(); });
	connect(&m_process, &QProcess::errorOccurred, this, [this](QProcess::ProcessError error) {
		if (error == QProcess::FailedToStart) m_timeout.stop();
	});
	connect(&m_process, &QProcess::finished, this, [this](int code, QProcess::ExitStatus status) {
		m_timeout.stop();
		if (status != QProcess::NormalExit || code != 0) return;
		m_version = availableVersion(m_process.readAllStandardOutput(), QCoreApplication::applicationVersion());
		emit changed();
	});
}

ApplicationUpdate::~ApplicationUpdate()
{
	if (m_process.state() != QProcess::NotRunning) {
		m_process.kill();
		m_process.waitForFinished();
	}
}

QString ApplicationUpdate::availableVersion(const QByteArray& output, const QString& currentVersion)
{
	const auto start = output.indexOf("<updates>");
	const auto end = output.indexOf("</updates>", start);
	if (start < 0 || end < 0) return {};
	QXmlStreamReader xml(output.mid(start, end + int(sizeof("</updates>") - 1) - start));
	QString version;
	while (!xml.atEnd()) {
		xml.readNext();
		if (xml.isStartElement() && xml.name() == u"update" && xml.attributes().value(u"id") == u"com.trajecta.app")
			version = xml.attributes().value(u"version").toString();
	}
	int suffix = 0;
	const auto candidate = QVersionNumber::fromString(version, &suffix);
	const auto current = QVersionNumber::fromString(currentVersion);
	if (xml.hasError() || suffix != version.size() || candidate.isNull() || current.isNull() ||
	    QVersionNumber::compare(candidate, current) <= 0)
		return {};
	return version;
}

QString ApplicationUpdate::toolPath() const
{
#ifdef Q_OS_WIN
	return QDir(QCoreApplication::applicationDirPath()).filePath(QStringLiteral("Uninstall Trajecta.exe"));
#elif defined(Q_OS_LINUX)
	return QDir(QCoreApplication::applicationDirPath()).filePath(QStringLiteral("../../TrajectaMaintenance"));
#else
	return {};
#endif
}

void ApplicationUpdate::check()
{
	const auto path = toolPath();
	if (path.isEmpty() || !QFileInfo::exists(path) || m_process.state() != QProcess::NotRunning) return;
	m_timeout.start(60000);
	m_process.start(path, {QStringLiteral("check-updates")});
}

bool ApplicationUpdate::openUpdater()
{
	if (!available()) return false;
	if (m_process.state() != QProcess::NotRunning) return false;
	if (QProcess::startDetached(toolPath(), {QStringLiteral("--start-updater")})) return true;
	m_error = tr("Could not start the update tool.");
	emit changed();
	return false;
}
