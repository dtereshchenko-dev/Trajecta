#include "SettingsStorage.h"

#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QMap>
#include <QSettings>
#include <QStandardPaths>
#include <QTemporaryDir>

namespace SettingsStorage {
namespace {
void removeEmptySource(QSettings& source)
{
	source.sync();
	const auto remaining = source.allKeys();
	if (source.status() != QSettings::NoError || !remaining.isEmpty())
		return;
	source.clear();
	source.sync();
	if (source.status() != QSettings::NoError)
		qWarning("Could not remove the empty previous settings store");
}
}

bool migrate(QSettings& source, const QString& destination, bool removeSource, bool* createdEmpty)
{
	if (createdEmpty)
		*createdEmpty = false;
	source.setFallbacksEnabled(false);
	const QFileInfo existing(destination);
	if (existing.exists()) {
		if (existing.isFile() && removeSource)
			removeEmptySource(source);
		return existing.isFile();
	}

	const QStringList keys = source.allKeys();
	if (source.status() != QSettings::NoError)
		return false;
	QMap<QString, QVariant> values;
	for (const QString& key : keys)
		values.insert(key, source.value(key));
	if (source.status() != QSettings::NoError)
		return false;

	const QFileInfo destinationInfo(destination);
	QDir directory = destinationInfo.absoluteDir();
	if (!directory.mkpath(QStringLiteral(".")))
		return false;

	QTemporaryDir temporary(directory.filePath(QStringLiteral(".settings-migration-XXXXXX")));
	if (!temporary.isValid())
		return false;
	const QString stagedPath = temporary.filePath(QStringLiteral("settings.ini"));
	QFile stagedFile(stagedPath);
	if (!stagedFile.open(QIODevice::WriteOnly))
		return false;
	stagedFile.close();
	{
		QSettings staged(stagedPath, QSettings::IniFormat);
		staged.setFallbacksEnabled(false);
		for (auto it = values.cbegin(); it != values.cend(); ++it)
			staged.setValue(it.key(), it.value());
		staged.sync();
		if (staged.status() != QSettings::NoError)
			return false;
	}
	{
		QSettings written(stagedPath, QSettings::IniFormat);
		written.setFallbacksEnabled(false);
		for (auto it = values.cbegin(); it != values.cend(); ++it) {
			QVariant restored = written.value(it.key());
			if (!written.contains(it.key()) ||
			    (restored.metaType() != it.value().metaType() && !restored.convert(it.value().metaType())) ||
			    restored != it.value()) {
				qWarning("Settings verification failed; the old settings were kept");
				return false;
			}
		}
		if (written.status() != QSettings::NoError)
			return false;
	}
	if (!QFile::rename(stagedPath, destination))
		return QFileInfo(destination).isFile();

	if (createdEmpty)
		*createdEmpty = keys.isEmpty();
	if (removeSource && !keys.isEmpty()) {
		for (auto it = values.cbegin(); it != values.cend(); ++it) {
			if (source.contains(it.key()) && source.value(it.key()) == it.value())
				source.remove(it.key());
		}
		source.sync();
		if (source.status() != QSettings::NoError)
			qWarning("Settings were copied to INI, but the old settings could not be removed");
	}
	if (removeSource)
		removeEmptySource(source);
	return true;
}

bool initialize(bool* createdEmpty)
{
	const QString organization = QCoreApplication::organizationName();
	const QString application = QCoreApplication::applicationName();
	QSettings source(QSettings::NativeFormat, QSettings::UserScope, organization, application);
#ifdef Q_OS_WIN
	const QString localData = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
	if (localData.isEmpty())
		return false;
	QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QFileInfo(localData).absolutePath());
#endif
	QSettings destination(QSettings::IniFormat, QSettings::UserScope, organization, application);
	const bool removeSource =
#ifdef Q_OS_WIN
	    true;
#else
	    false;
#endif
	if (!migrate(source, destination.fileName(), removeSource, createdEmpty))
		return false;
	QSettings::setDefaultFormat(QSettings::IniFormat);
	return true;
}
}
