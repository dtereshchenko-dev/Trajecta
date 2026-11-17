#pragma once

#include <QFileInfo>
#include <QSettings>
#include <QScopeGuard>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QTest>
#include <QUuid>
#include "settings/FitSettingsAggregator.h"
#include "settings/SettingsStorage.h"

#ifdef Q_OS_WIN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

class TestSettingsStorage : public QObject {
	Q_OBJECT

   private slots:

	void migratesAllValuesOnce()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings source(directory.filePath("old.ini"), QSettings::IniFormat);
		source.setValue("optimizer/integratorInitialStepDays", 2.5);
		source.setValue("appearance/theme", "slate");
		source.setValue("recentFiles", QStringList{"one.xml", "two.psv"});
		source.sync();
		const QString path = directory.filePath("config/Trajecta.ini");
		QVERIFY(SettingsStorage::migrate(source, path));
		QVERIFY(QFileInfo::exists(path));
		QSettings result(path, QSettings::IniFormat);
		QCOMPARE(result.value("optimizer/integratorInitialStepDays").toDouble(), 2.5);
		QCOMPARE(result.value("appearance/theme").toString(), QString("slate"));
		QCOMPARE(result.value("recentFiles").toStringList(), (QStringList{"one.xml", "two.psv"}));

		result.setValue("appearance/theme", "graphite");
		result.sync();
		source.setValue("appearance/theme", "light");
		source.sync();
		QVERIFY(SettingsStorage::migrate(source, path));
		QSettings reopened(path, QSettings::IniFormat);
		QCOMPARE(reopened.value("appearance/theme").toString(), QString("graphite"));
		QCOMPARE(source.value("appearance/theme").toString(), QString("light"));
	}

	void createsEmptySettingsForNewUsers()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings source(directory.filePath("missing.ini"), QSettings::IniFormat);
		const QString path = directory.filePath("new.ini");
		bool createdEmpty = false;
		QVERIFY(SettingsStorage::migrate(source, path, true, &createdEmpty));
		QVERIFY(createdEmpty);
		QVERIFY(QFileInfo::exists(path));
		QSettings result(path, QSettings::IniFormat);
		QVERIFY(result.allKeys().isEmpty());
	}

	void writesFitDefaultsOnFirstRun()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const auto previousFormat = QSettings::defaultFormat();
		const QSettings previous(QSettings::IniFormat, QSettings::UserScope, QCoreApplication::organizationName(),
		                         QCoreApplication::applicationName());
		const QString originalRoot = QFileInfo(QFileInfo(previous.fileName()).absolutePath()).absolutePath();
		const auto restore = qScopeGuard([&] {
			QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, originalRoot);
			QSettings::setDefaultFormat(previousFormat);
		});
		QSettings::setDefaultFormat(QSettings::IniFormat);
		QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
		FitSettingsAggregator settings;
		settings.load();
		settings.save();
		QSettings saved;
		QVERIFY(QFileInfo::exists(saved.fileName()));
		QVERIFY(saved.value("integrator/nBodyPropagation").toBool());
		QCOMPARE(saved.value("integrator/initialStepDays").toDouble(), 10.0);
	}

	void removesOldValuesAfterSuccessfulCopy()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings source(directory.filePath("old.ini"), QSettings::IniFormat);
		source.setValue("optimizer/pertEnabled", true);
		source.setValue("optimizer/integratorInitialStepDays", 2.5);
		source.setValue("recentFiles", QStringList{"one.xml", "two.psv"});
		source.sync();
		const QString path = directory.filePath("new.ini");
		bool createdEmpty = true;
		QVERIFY(SettingsStorage::migrate(source, path, true, &createdEmpty));
		QVERIFY(!createdEmpty);
		QVERIFY(source.allKeys().isEmpty());
		QSettings result(path, QSettings::IniFormat);
		QVERIFY(result.value("optimizer/pertEnabled").toBool());
		QCOMPARE(result.value("optimizer/integratorInitialStepDays").toDouble(), 2.5);
		QCOMPARE(result.value("recentFiles").toStringList(), (QStringList{"one.xml", "two.psv"}));
	}

	void keepsOldValuesWhenDestinationAlreadyExists()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings source(directory.filePath("old.ini"), QSettings::IniFormat);
		source.setValue("appearance/theme", "slate");
		source.sync();
		const QString path = directory.filePath("new.ini");
		QSettings result(path, QSettings::IniFormat);
		result.setValue("appearance/theme", "graphite");
		result.sync();
		QVERIFY(SettingsStorage::migrate(source, path, true));
		QCOMPARE(source.value("appearance/theme").toString(), QString("slate"));
		QCOMPARE(result.value("appearance/theme").toString(), QString("graphite"));
	}

#ifdef Q_OS_WIN
	void windowsIniPathIsInsideApplicationDataDirectory()
	{
		const QString originalOrganization = QCoreApplication::organizationName();
		const QString originalApplication = QCoreApplication::applicationName();
		const QSettings previous(QSettings::IniFormat, QSettings::UserScope, originalOrganization, originalApplication);
		const QString originalRoot = QFileInfo(QFileInfo(previous.fileName()).absolutePath()).absolutePath();
		const auto restore = qScopeGuard([&] {
			QCoreApplication::setOrganizationName(originalOrganization);
			QCoreApplication::setApplicationName(originalApplication);
			QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, originalRoot);
		});
		QCoreApplication::setOrganizationName(QStringLiteral("Trajecta"));
		QCoreApplication::setApplicationName(QStringLiteral("Trajecta"));
		const QString localData = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
		QVERIFY(!localData.isEmpty());
		QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QFileInfo(localData).absolutePath());
		const QSettings settings(QSettings::IniFormat, QSettings::UserScope, QStringLiteral("Trajecta"),
		                         QStringLiteral("Trajecta"));
		QCOMPARE(QDir::cleanPath(settings.fileName()), QDir::cleanPath(QDir(localData).filePath("Trajecta.ini")));
	}

	void migratesIsolatedWindowsRegistryValues()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const QString name =
		    QStringLiteral("SettingsMigrationTest-") + QUuid::createUuid().toString(QUuid::WithoutBraces);
		QSettings source(QSettings::NativeFormat, QSettings::UserScope, name, QStringLiteral("Test"));
		source.setFallbacksEnabled(false);
		const auto cleanup = qScopeGuard([&source] {
			source.clear();
			source.sync();
		});
		source.setValue("optimizer/pertEnabled", true);
		source.setValue("optimizer/planetMask", 6399u);
		source.setValue("optimizer/integratorInitialStepDays", 2.5);
		source.setValue("recentFiles", QStringList{"one.xml", "two.psv"});
		source.sync();
		QCOMPARE(source.status(), QSettings::NoError);
		const QString path = directory.filePath("new.ini");
		QVERIFY(SettingsStorage::migrate(source, path, true));
		QVERIFY(source.allKeys().isEmpty());
		QSettings result(path, QSettings::IniFormat);
		QVERIFY(result.value("optimizer/pertEnabled").toBool());
		QCOMPARE(result.value("optimizer/planetMask").toUInt(), 6399u);
		QCOMPARE(result.value("optimizer/integratorInitialStepDays").toDouble(), 2.5);
		QCOMPARE(result.value("recentFiles").toStringList(), (QStringList{"one.xml", "two.psv"}));
	}

	void removesEmptyRegistrySubkeysWhenIniAlreadyExists()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const QString name =
		    QStringLiteral("SettingsMigrationTest-") + QUuid::createUuid().toString(QUuid::WithoutBraces);
		QSettings source(QSettings::NativeFormat, QSettings::UserScope, name, QStringLiteral("Test"));
		source.setFallbacksEnabled(false);
		const auto cleanup = qScopeGuard([&source] {
			source.clear();
			source.sync();
		});
		const QString optimizerKey = QStringLiteral("Software\\%1\\Test\\optimizer").arg(name);
		const QString childKey = optimizerKey + QStringLiteral("\\photometry");
		HKEY key = nullptr;
		QCOMPARE(RegCreateKeyExW(HKEY_CURRENT_USER, reinterpret_cast<LPCWSTR>(childKey.utf16()), 0, nullptr, 0,
		                         KEY_WRITE, nullptr, &key, nullptr),
		         LSTATUS(ERROR_SUCCESS));
		RegCloseKey(key);
		const QString path = directory.filePath("new.ini");
		QSettings result(path, QSettings::IniFormat);
		result.setValue("fit/useWeights", true);
		result.sync();
		QVERIFY(SettingsStorage::migrate(source, path, true));
		key = nullptr;
		const LSTATUS status =
		    RegOpenKeyExW(HKEY_CURRENT_USER, reinterpret_cast<LPCWSTR>(optimizerKey.utf16()), 0, KEY_READ, &key);
		if (status == ERROR_SUCCESS)
			RegCloseKey(key);
		QCOMPARE(status, LSTATUS(ERROR_FILE_NOT_FOUND));
		QVERIFY(result.value("fit/useWeights").toBool());
	}
#endif

	void rejectsDirectoryAsDestination()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings source(directory.filePath("old.ini"), QSettings::IniFormat);
		source.setValue("appearance/theme", "slate");
		source.sync();
		QVERIFY(!SettingsStorage::migrate(source, directory.path(), true));
		QCOMPARE(source.value("appearance/theme").toString(), QString("slate"));
	}
};
