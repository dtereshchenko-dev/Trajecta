#include <QTest>
#include <QApplication>
#include <QProcess>
#include <QProcessEnvironment>
#include <QSettings>
#include <QTemporaryDir>

#include "tst_recentfiles.h"
#include "tst_orbitgeometry.h"
#include "tst_cometorbits.h"
#include "tst_gridgeometry.h"
#include "tst_diagnosticsmodel.h"
#include "tst_nongravsettings.h"
#include "tst_orbitengine.h"
#include "tst_workspace.h"

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	app.setOrganizationName("TestOrg");
	app.setApplicationName("TestApp");
	if (qEnvironmentVariableIsSet("ORBITAL_DIAGNOSTICS_ONLY")) {
		TestDiagnosticsModel test;
		return QTest::qExec(&test, argc, argv);
	}

	int status = 0;
	if (qEnvironmentVariableIsSet("ORBITAL_WORKSPACE_ONLY")) {
		TestWorkspace test;
		return QTest::qExec(&test, argc, argv);
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
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestOrbitGeometry test;
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestGridGeometry test;
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestOrbitEngine test;
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestDiagnosticsModel test;
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestNonGravSettings test;
		status |= QTest::qExec(&test, argc, argv);
	}
	{
		TestCometOrbits test;
		status |= QTest::qExec(&test, argc, argv);
	}

	return status;
}
