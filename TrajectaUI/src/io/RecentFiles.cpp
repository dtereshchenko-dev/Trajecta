#include "io/RecentFiles.h"
#include <algorithm>
#include <QSettings>
#include <QFileInfo>
#include <QDir>

RecentFiles::RecentFiles(QObject* parent) : QObject(parent) { load(); }

QStringList RecentFiles::files() const { return m_files; }

int RecentFiles::count() const { return m_files.size(); }

void RecentFiles::addFile(const QUrl& url)
{
	QString path = url.toLocalFile();
	if (path.isEmpty())
		return;

	m_files.removeAll(path);
	m_files.prepend(path);

	while (m_files.size() > MaxRecent)
		m_files.removeLast();

	save();
	emit filesChanged();
}

void RecentFiles::removeFile(const QString& path)
{
	if (m_files.removeAll(path) > 0) {
		save();
		emit filesChanged();
	}
}

void RecentFiles::clear()
{
	m_files.clear();
	save();
	emit filesChanged();
}

QString RecentFiles::fileName(int index) const
{
	if (index < 0 || index >= m_files.size())
		return {};

	return QFileInfo(m_files[index]).fileName();
}

QString RecentFiles::filePath(int index) const
{
	if (index < 0 || index >= m_files.size())
		return {};

	return m_files[index];
}

QString RecentFiles::displayPath(int index) const
{
	if (index < 0 || index >= m_files.size())
		return {};

	return QDir::toNativeSeparators(m_files[index]);
}

void RecentFiles::load()
{
	QSettings settings;
	m_files = settings.value("recentFiles").toStringList();

	m_files.erase(
	    std::remove_if(m_files.begin(), m_files.end(), [](const QString& f) { return !QFileInfo::exists(f); }),
	    m_files.end());
}

void RecentFiles::save()
{
	QSettings settings;
	settings.setValue("recentFiles", m_files);
}
