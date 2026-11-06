#include <QApplication>
#include <QQmlApplicationEngine>
#include <QQuickStyle>
#include <QIcon>
#include <QStandardPaths>
#include <QDir>
#include "settings/FitSettingsAggregator.h"
#include "orbital_core/Types.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/RuntimeData.h"
#include "logger/Logger.h"
#include "app/ApplicationSession.h"
#include "features/settings/Sb441Download.h"
#include <QFileInfo>
#include <memory>
#ifdef Q_OS_WIN
#include "app/TaskbarProgress.h"
#include <QWindow>
#endif

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	app.setOrganizationName("Trajecta");
	app.setApplicationName("Trajecta");
	app.setApplicationDisplayName("Trajecta");
	app.setApplicationVersion(QStringLiteral(TRAJECTA_VERSION));

	getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();

	auto logDir = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
	if (!logDir.isEmpty()) {
		orbital_core::setObservatoryDataDirectory(QDir(logDir).filePath("observatories").toStdWString());
		orbital_core::setEphemerisDataDirectory(QDir(logDir).filePath("ephemerides").toStdWString());
	}
	QDir().mkpath(logDir + "/log");
	Logger::setLogFile((logDir + "/log/process.log").toStdString());

	DE440Source de440;
	setDefaultEphemerisSource(&de440);

	std::unique_ptr<SpkReader> sb441;
	const auto downloadedSb441 = Sb441Download::installedPath();
	if (QFileInfo::exists(downloadedSb441))
		sb441 = std::make_unique<SpkReader>(downloadedSb441.toUtf8().toStdString());
	if (!sb441 || !sb441->isAvailable())
		sb441 = std::make_unique<SpkReader>();
	setDefaultAsteroidSource(sb441.get());

	QQuickStyle::setStyle("Basic");
	app.setWindowIcon(QIcon(":/qt/qml/TrajectaUI/resources/icons/workspace/app-logo.svg"));

	FitSettingsAggregator settings;
	settings.load();

	ApplicationSession session(settings);
	session.applicationUpdate()->check();
	session.dataSources()->observatoryUpdate()->check();
	QQmlApplicationEngine engine;
	engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});

	QObject::connect(
	    &engine, &QQmlApplicationEngine::objectCreationFailed, &app, []() { QCoreApplication::exit(-1); },
	    Qt::QueuedConnection);
	engine.loadFromModule("TrajectaUI", "WorkspaceWindow");

#ifdef Q_OS_WIN
	std::unique_ptr<TaskbarProgress> taskbar;
	QMetaObject::Connection taskbarConnection;
	if (!engine.rootObjects().isEmpty()) {
		if (auto* window = qobject_cast<QWindow*>(engine.rootObjects().front())) {
			taskbar = std::make_unique<TaskbarProgress>(window);
			taskbarConnection = QObject::connect(&session, &ApplicationSession::changed, &app,
			                                     [&session, &taskbar] { taskbar->setBusy(session.busy()); });
		}
	}
#endif

	const int result = app.exec();
#ifdef Q_OS_WIN
	QObject::disconnect(taskbarConnection);
#endif
	return result;
}
