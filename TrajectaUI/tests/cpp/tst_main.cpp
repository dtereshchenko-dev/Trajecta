#include <QTest>
#include <QApplication>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QProcessEnvironment>
#include <QSettings>
#include <QTemporaryDir>
#include <cstdio>
#include "orbital_core/RuntimeData.h"

#include "tst_recentfiles.h"
#include "tst_orbitgeometry.h"
#include "tst_cometorbits.h"
#include "tst_gridgeometry.h"
#include "tst_diagnosticsmodel.h"
#include "tst_nongravsettings.h"
#include "tst_settingsstorage.h"
#include "tst_orbitengine.h"
#include "tst_workspace.h"

static int runTest(QObject* test, int argc, char* argv[])
{
	const QString logDirectory = qEnvironmentVariable("TRAJECTA_QT_TEST_LOG_DIR");
	if (logDirectory.isEmpty())
		return QTest::qExec(test, argc, argv);

	if (!QDir().mkpath(logDirectory)) {
		std::fprintf(stderr, "Cannot create Qt test log directory: %s\n", qPrintable(logDirectory));
		return 1;
	}
	const QString logPath = QDir(logDirectory)
	                            .filePath(QStringLiteral("%1-%2.txt")
	                                          .arg(QString::fromLatin1(test->metaObject()->className()))
	                                          .arg(QCoreApplication::applicationPid()));
	QStringList arguments = QCoreApplication::arguments();
	arguments << QStringLiteral("-o") << logPath + QStringLiteral(",txt");
	const int status = QTest::qExec(test, arguments);
	if (status != 0) {
		QFile log(logPath);
		if (log.open(QIODevice::ReadOnly)) {
			const QByteArray output = log.readAll();
			std::fwrite(output.constData(), 1, static_cast<size_t>(output.size()), stderr);
			std::fflush(stderr);
		} else {
			std::fprintf(stderr, "Cannot read Qt test log: %s\n", qPrintable(logPath));
		}
	}
	return status;
}

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	app.setOrganizationName("TestOrg");
	app.setApplicationName("TestApp");
	orbital_core::setObservatoryDataDirectory(
	    QDir(QStringLiteral(ORBITAL_UI_SOURCE_DIR)).filePath("../data").toStdWString());
	orbital_core::setEphemerisDataDirectory(
	    QDir(QStringLiteral(ORBITAL_UI_SOURCE_DIR)).filePath("../data").toStdWString());
	orbital_core::setEarthPckDataDirectory(
	    QDir(QStringLiteral(ORBITAL_UI_SOURCE_DIR)).filePath("../data/earth_orientation").toStdWString());
	if (qEnvironmentVariableIsSet("ORBITAL_DIAGNOSTICS_ONLY")) {
		TestDiagnosticsModel test;
		return runTest(&test, argc, argv);
	}
	if (qEnvironmentVariableIsSet("ORBITAL_SETTINGS_STORAGE_ONLY")) {
		TestSettingsStorage test;
		return runTest(&test, argc, argv);
	}

	int status = 0;
	if (qEnvironmentVariableIsSet("ORBITAL_WORKSPACE_ONLY")) {
		TestWorkspace test;
		return runTest(&test, argc, argv);
	}
	if (!qEnvironmentVariableIsSet("ORBITAL_MODELS_ONLY")) {
		QProcess workspace;
		auto environment = QProcessEnvironment::systemEnvironment();
		environment.insert("ORBITAL_WORKSPACE_ONLY", "1");
		environment.insert("QT_QUICK_CONTROLS_STYLE", "Basic");
		workspace.setProcessEnvironment(environment);
		workspace.setProcessChannelMode(QProcess::ForwardedChannels);
		workspace.start(QCoreApplication::applicationFilePath(), QCoreApplication::arguments().mid(1));
		if (!workspace.waitForFinished(180000)) {
			workspace.kill();
			workspace.waitForFinished();
			qCritical("Workspace tests did not finish: %s", qPrintable(workspace.errorString()));
			status |= 1;
		} else {
			status |= workspace.exitStatus() == QProcess::NormalExit ? workspace.exitCode() : 1;
		}
	}

	QTemporaryDir settingsDirectory;
	if (!settingsDirectory.isValid())
		return 1;
	QSettings::setDefaultFormat(QSettings::IniFormat);
	QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, settingsDirectory.path());

	{
		TestRecentFiles test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestOrbitGeometry test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestGridGeometry test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestOrbitEngine test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestDiagnosticsModel test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestNonGravSettings test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestSettingsStorage test;
		status |= runTest(&test, argc, argv);
	}
	{
		TestCometOrbits test;
		status |= runTest(&test, argc, argv);
	}

	return status;
}
