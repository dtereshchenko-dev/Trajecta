#pragma once

#include <QTest>
#include <QElapsedTimer>
#include <QGraphicsScene>
#include <array>
#include <QAbstractItemModelTester>
#include <QTemporaryDir>
#include <QQmlApplicationEngine>
#include <QQmlComponent>
#include <QQmlListReference>
#include <QQuickWindow>
#include <QQuickItem>
#include <QQuickItemGrabResult>
#include <QDirIterator>
#include <QLibraryInfo>
#include <QScopeGuard>
#include <QMimeData>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QIcon>
#include <QImageReader>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTcpServer>
#include <QTcpSocket>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include "app/ApplicationSession.h"
#include "features/solution/ResidualMarker.h"
#include "geometry/OrbitGeometry.h"
#include "geometry/GridGeometry.h"
#include "util/FormatUtils.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include <orbital_core/CometaryCoordinates.h>
#include <orbital_core/TimeScale.h>
#include <orbital_core/hazard/EquinoctialCoordinates.h>
#include <orbital_core/hazard/BodyFixedFrame.h>

class TestWorkspace : public QObject {
	Q_OBJECT

	QTemporaryDir m_settings;

	static QQuickItem* findItem(QQuickItem* parent, const QString& name)
	{
		if (parent->objectName() == name)
			return parent;
		for (auto* child : parent->childItems())
			if (auto* found = findItem(child, name))
				return found;
		return nullptr;
	}

	static orbital_core::FitPipelineResult fixture()
	{
		orbital_core::FitPipelineResult result;
		result.observations.objectName = "Fixture";
		orbital_core::hazard::equinoctial::Vector6 elements;
		elements << 1.4, 0.12, 0.16, 0.2, -0.14, 1.3;
		const auto state = orbital_core::hazard::equinoctial::toCartesian(elements).value();
		result.fit.orbit = make_orbit_from_states_vectors(
		    StateVectors(Vec3D(state(0), state(1), state(2)), Vec3D(state(3), state(4), state(5)) / Constants::GAUSS_K),
		    2460000.5);
		result.fit.fitEpochOrbit = std::make_shared<Orbit>(*result.fit.orbit);
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(6, 6) * 1e-8;
		result.fit.converged = true;
		result.fit.cost = 0.2;
		for (int i = 0; i < 3; ++i) {
			ObsNode observation;
			observation.JD = 2460000.5 + i;
			observation.RA = 1;
			observation.DEC = -0.005;
			observation.obsCode = i == 1 ? "B02" : "A01";
			result.observations.all.push_back(observation);
			Residual residual{};
			residual.dRA_arcsec = i * 0.1;
			residual.dDEC_arcsec = -i * 0.1;
			result.fit.allResiduals.add(residual);
			if (i != 1)
				result.fit.residuals.add(residual);
		}
		result.fit.keptObsIndices = {0, 2};
		return result;
	}

	static void registerQml()
	{
		const auto root = QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/";
		qmlRegisterUncreatableType<ApplicationSession>("TrajectaUI", 1, 0, "ApplicationSession", "Application-owned");
		qmlRegisterUncreatableType<CalculationSettingsModel>("TrajectaUI", 1, 0, "CalculationSettingsModel",
		                                                     "Session-owned");
		qmlRegisterUncreatableType<DataSourcesModel>("TrajectaUI", 1, 0, "DataSourcesModel", "Session-owned");
		qmlRegisterUncreatableType<Sb441Download>("TrajectaUI", 1, 0, "Sb441Download", "Session-owned");
		qmlRegisterUncreatableType<ObservatoryUpdate>("TrajectaUI", 1, 0, "ObservatoryUpdate", "Session-owned");
		qmlRegisterUncreatableType<EarthPckUpdate>("TrajectaUI", 1, 0, "EarthPckUpdate", "Session-owned");
		qmlRegisterUncreatableType<ApplicationUpdate>("TrajectaUI", 1, 0, "ApplicationUpdate", "Session-owned");
		qmlRegisterUncreatableType<SolutionViewModel>("TrajectaUI", 1, 0, "SolutionViewModel", "Session-owned");
		qmlRegisterUncreatableType<ObservationTableModel>("TrajectaUI", 1, 0, "ObservationTableModel",
		                                                  "Solution-owned");
		qmlRegisterUncreatableType<StationResidualModel>("TrajectaUI", 1, 0, "StationResidualModel",
		                                                 "Observation-owned");
		qmlRegisterUncreatableType<CovarianceTableModel>("TrajectaUI", 1, 0, "CovarianceTableModel", "Solution-owned");
		qmlRegisterUncreatableType<DiagnosticsModel>("TrajectaUI", 1, 0, "DiagnosticsModel", "Solution-owned");
		qmlRegisterUncreatableType<EphemerisViewModel>("TrajectaUI", 1, 0, "EphemerisViewModel", "Session-owned");
		qmlRegisterUncreatableType<CloseApproachesViewModel>("TrajectaUI", 1, 0, "CloseApproachesViewModel",
		                                                     "Session-owned");
		qmlRegisterUncreatableType<CloseApproachTableModel>("TrajectaUI", 1, 0, "CloseApproachTableModel",
		                                                    "Approaches-owned");
		qmlRegisterUncreatableType<EphemerisTableModel>("TrajectaUI", 1, 0, "EphemerisTableModel", "Ephemeris-owned");
		qmlRegisterUncreatableType<RiskMonitorViewModel>("TrajectaUI", 1, 0, "RiskMonitorViewModel", "Session-owned");
		qmlRegisterUncreatableType<MovViewModel>("TrajectaUI", 1, 0, "MovViewModel", "Session-owned");
		qmlRegisterUncreatableType<AtmosphereViewModel>("TrajectaUI", 1, 0, "AtmosphereViewModel", "Session-owned");
		qmlRegisterUncreatableType<RiskEventTableModel>("TrajectaUI", 1, 0, "RiskEventTableModel", "Monitor-owned");
		qmlRegisterType<ResidualMarker>("TrajectaUI", 1, 0, "ResidualMarker");
		qmlRegisterUncreatableType<OrbitModel>("TrajectaUI", 1, 0, "OrbitModel", "Solution-owned");
		qmlRegisterUncreatableType<OrbitSceneModel>("TrajectaUI", 1, 0, "OrbitSceneModel", "Solution-owned");
		qmlRegisterType<OrbitGeometry>("TrajectaUI", 1, 0, "OrbitGeometry");
		qmlRegisterType<GridGeometry>("TrajectaUI", 1, 0, "GridGeometry");
		qmlRegisterUncreatableType<RecentFiles>("TrajectaUI", 1, 0, "RecentFiles", "Session-owned");
		qmlRegisterSingletonType(QUrl::fromLocalFile(root + "shared/AppTheme.qml"), "TrajectaUI", 1, 0, "AppTheme");
		qmlRegisterSingletonType(QUrl::fromLocalFile(root + "shared/DialogState.qml"), "TrajectaUI", 1, 0,
		                         "DialogState");
		qmlRegisterSingletonType(QUrl::fromLocalFile(root + "shared/BodyPalette.qml"), "TrajectaUI", 1, 0,
		                         "BodyPalette");
		for (const QString directory :
		     {"controls", "shell", "components/shared", "features/solution", "features/observations",
		      "features/settings", "features/orbit", "features/ephemerides", "features/risk", "features/atmosphere",
		      "features/encounters"}) {
			QDirIterator files(root + directory, {"*.qml"}, QDir::Files);
			while (files.hasNext()) {
				const auto path = files.next();
				const auto name = QFileInfo(path).baseName().toUtf8();
				qmlRegisterType(QUrl::fromLocalFile(path), "TrajectaUI", 1, 0, name.constData());
			}
		}
	}

	static QJsonObject atmosphericArchiveFixture()
	{
		namespace hz = orbital_core::hazard;
		const double tt = 60000.25;
		const double tdb = orbital_core::time::toMjdTdb(tt, orbital_core::time::Scale::TT);
		const auto r = strewn::ecefPosition({0.5, 0.3, 100000});
		const auto velocity = strewn::enuToEcef({0, 500, -1000}, 0.5, 0.3);
		const auto frame = hz::bodyFixedTransformMjd(hz::Planet::Earth, tdb);
		const auto inverse = [&frame](const Vec3D& v) {
			return Helpers::ecvToEcl(frame.equatorialRows[0] * v.x + frame.equatorialRows[1] * v.y +
			                         frame.equatorialRows[2] * v.z);
		};
		const auto inertial = inverse({r[0] / 1000, r[1] / 1000, r[2] / 1000});
		const auto rate = hz::bodyFixedVelocityMjd(hz::Planet::Earth, inertial, Vec3D{}, tdb);
		const auto v = inverse(Vec3D{velocity[0] / 1000, velocity[1] / 1000, velocity[2] / 1000} - rate);
		const QJsonObject entry{{"epoch", tt},
		                        {"relative_position_km", QJsonArray{inertial.x, inertial.y, inertial.z}},
		                        {"relative_velocity_km_s", QJsonArray{v.x, v.y, v.z}},
		                        {"latitude_deg", 0.5 * 180 / Constants::PI},
		                        {"longitude_deg", 0.3 * 180 / Constants::PI}};
		return {{"schema", "orbital-interval-mc-3"},
		        {"time_scale", "MJD TT"},
		        {"frame", "heliocentric ecliptic J2000"},
		        {"state_units", "AU, AU/day"},
		        {"target", 2},
		        {"altitude_km", 100},
		        {"body_fixed_frame", QString::fromStdString(std::string(frame.frame))},
		        {"body_fixed_model", QString::fromStdString(hz::earthBodyFixedFrameModel())},
		        {"method", "covariance"},
		        {"requested", 2},
		        {"completed", 2},
		        {"cancelled", false},
		        {"failed_fit", 0},
		        {"failed_propagation", 0},
		        {"cancelled_clones", 0},
		        {"unstarted", 0},
		        {"entry_clones", 1},
		        {"start", tt - 1},
		        {"end", tt + 1},
		        {"clones",
		         QJsonArray{
		             QJsonObject{{"index", 0}, {"status", 1}, {"entry_before_interval", false}, {"entry", entry}},
		             QJsonObject{
		                 {"index", 1}, {"status", 1}, {"entry_before_interval", false}, {"entry", QJsonValue::Null}}}},
		        {"events", QJsonArray{QJsonObject{{"start", tt}, {"end", tt}, {"entries", 1}}}}};
	}

   private slots:

	void applicationUpdateParsesMaintenanceToolOutput()
	{
		const QByteArray output =
		    "[12] Loading component scripts...\n<?xml version=\"1.0\"?>\n"
		    "<updates><update id=\"other\" version=\"9.0.0\"/>"
		    "<update id=\"com.trajecta.app\" version=\"1.0.2\"/></updates>\n";
		QCOMPARE(ApplicationUpdate::availableVersion(output, "1.0.1"), QStringLiteral("1.0.2"));
		QVERIFY(ApplicationUpdate::availableVersion(output, "1.0.2").isEmpty());
		QVERIFY(ApplicationUpdate::availableVersion(output, "2.0.0").isEmpty());
		QVERIFY(ApplicationUpdate::availableVersion("network error", "1.0.1").isEmpty());
		QVERIFY(ApplicationUpdate::availableVersion(
		            "<updates><update id=\"com.trajecta.app\" version=\"1.0.2\"></updates>", "1.0.1")
		            .isEmpty());
		QVERIFY(ApplicationUpdate::availableVersion(
		            "<updates><update id=\"com.trajecta.app\" version=\"1.0.2-beta\"/></updates>", "1.0.1")
		            .isEmpty());
	}

	void sb441DownloadVerifiesBeforeReplacing()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const QByteArray payload("local test payload");
		QTcpServer server;
		QVERIFY(server.listen(QHostAddress::LocalHost));
		connect(&server, &QTcpServer::newConnection, &server, [&] {
			auto* socket = server.nextPendingConnection();
			connect(socket, &QTcpSocket::readyRead, socket, [socket, payload] {
				auto request = socket->property("request").toByteArray() + socket->readAll();
				socket->setProperty("request", request);
				if (!request.contains("\r\n\r\n") || socket->property("responded").toBool())
					return;
				socket->setProperty("responded", true);
				socket->write("HTTP/1.1 200 OK\r\nContent-Length: " + QByteArray::number(payload.size()) +
				              "\r\nConnection: close\r\n\r\n" + payload);
				socket->disconnectFromHost();
			});
		});
		const auto url = QUrl(QString("http://127.0.0.1:%1/kernel").arg(server.serverPort()));
		const auto target = directory.filePath("sb441-n16.bsp");
		QFile old(target);
		QVERIFY(old.open(QIODevice::WriteOnly));
		QCOMPARE(old.write("original"), 8);
		old.close();
		Sb441Download invalid(url, directory.path(), payload.size(), QByteArray(32, '\0'));
		QSignalSpy rejected(&invalid, &Sb441Download::installed);
		invalid.start();
		QTRY_VERIFY_WITH_TIMEOUT(!invalid.running(), 10000);
		QCOMPARE(rejected.count(), 0);
		QVERIFY(invalid.message().contains("verification", Qt::CaseInsensitive));
		QVERIFY(old.open(QIODevice::ReadOnly));
		QCOMPARE(old.readAll(), QByteArray("original"));
		old.close();
		const auto hash = QCryptographicHash::hash(payload, QCryptographicHash::Sha256);
		Sb441Download valid(url, directory.path(), payload.size(), hash);
		QSignalSpy installed(&valid, &Sb441Download::installed);
		valid.start();
		QTRY_COMPARE_WITH_TIMEOUT(installed.count(), 1, 10000);
		QVERIFY(old.open(QIODevice::ReadOnly));
		QCOMPARE(old.readAll(), payload);
	}

	void sb441DownloadCancellationPreservesExistingFile()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QFile old(directory.filePath("sb441-n16.bsp"));
		QVERIFY(old.open(QIODevice::WriteOnly));
		QCOMPARE(old.write("original"), 8);
		old.close();
		QTcpServer server;
		QVERIFY(server.listen(QHostAddress::LocalHost));
		connect(&server, &QTcpServer::newConnection, &server, [&] {
			auto* socket = server.nextPendingConnection();
			connect(socket, &QTcpSocket::readyRead, socket, [socket] {
				auto request = socket->property("request").toByteArray() + socket->readAll();
				socket->setProperty("request", request);
				if (!request.contains("\r\n\r\n") || socket->property("responded").toBool())
					return;
				socket->setProperty("responded", true);
				socket->write("HTTP/1.1 200 OK\r\nContent-Length: 10000\r\n\r\npartial");
				socket->flush();
			});
		});
		const auto url = QUrl(QString("http://127.0.0.1:%1/kernel").arg(server.serverPort()));
		Sb441Download download(url, directory.path(), 10000, QByteArray(32, '\0'));
		QSignalSpy installed(&download, &Sb441Download::installed);
		download.start();
		QTRY_VERIFY_WITH_TIMEOUT(download.progress() > 0, 10000);
		download.cancel();
		QTRY_VERIFY_WITH_TIMEOUT(!download.running(), 10000);
		QCOMPARE(installed.count(), 0);
		QVERIFY(old.open(QIODevice::ReadOnly));
		QCOMPARE(old.readAll(), QByteArray("original"));
	}

	void solutionStatusDistinguishesSelectionFromCorrection()
	{
		using Status = orbital_core::FitSelectionStatus;
		SolutionViewModel solution;
		auto result = fixture();
		result.fit.converged = false;
		for (const auto& [status, text] :
		     std::array{std::pair{Status::InsufficientObservations, "Too few observations after rejection"},
		                std::pair{Status::PassLimit, "Observation selection limit reached"},
		                std::pair{Status::SafeguardRollback, "Observation selection safeguard reached"},
		                std::pair{Status::CorrectionFailed, "Fit not converged"}}) {
			result.fit.selectionStatus = status;
			QVERIFY(solution.publish(result, {}, {}));
			QCOMPARE(solution.status(), QString::fromLatin1(text));
		}
		result.fit.converged = true;
		result.fit.selectionStatus = Status::Stable;
		QVERIFY(solution.publish(result, {}, {}));
		QCOMPARE(solution.status(), QStringLiteral("Converged"));
	}

	void closeApproachTextColumnsAlign()
	{
		namespace hz = orbital_core::hazard;
		hz::NominalCloseApproach earth, moon, jupiter;
		earth.mjdTdb = 60000.25;
		earth.encounter.planet = hz::Planet::Earth;
		earth.encounter.distAu = 0.000415351382;
		earth.minimumDistanceAu = 1.05756781e-5;
		earth.maximumDistanceAu = 0.00141083247;
		earth.timeUncertaintyDays = 0.5 / 1440;
		earth.relativeSpeedKmS = 6.24540961;
		earth.encounter.vInfKmS = 5.11616459;
		moon = earth;
		moon.mjdTdb += 1;
		moon.encounter.planet = hz::Planet::Moon;
		moon.minimumDistanceAu.reset();
		moon.timeUncertaintyDays = 123456.25;
		jupiter = earth;
		jupiter.mjdTdb += 2;
		jupiter.encounter.planet = hz::Planet::Jupiter;
		jupiter.timeUncertaintyDays.reset();
		jupiter.encounter.vInfKmS = 0;
		CloseApproachTableModel table;
		table.setRows({earth, moon, jupiter});
		for (int mode = 0; mode < 3; ++mode) {
			if (mode == 1) {
				table.setDistanceUnit(1);
				table.setVelocityUnit(1);
				table.toggleSort(0);
			}
			if (mode == 2)
				table.setRows({});
			const auto text = table.toPlainText();
			QVERIFY(!text.contains('\t'));
			const auto lines = text.split('\n');
			QCOMPARE(lines.size(), table.count() + 2);
			QVERIFY(lines.last().isEmpty());
			for (int row = -1; row < table.count(); ++row) {
				const auto& line = lines[row + 1];
				QCOMPARE(line.size(), lines.first().size());
				qsizetype offset = 0;
				for (int c = 0; c < table.columnCount(); ++c) {
					const auto heading = table.header(c).replace('\n', ' ');
					auto width = heading.size();
					for (int r = 0; r < table.count(); ++r)
						width = std::max(width, table.cell(r, c).size());
					const auto value = row < 0 ? heading : table.cell(row, c);
					QCOMPARE(line.mid(offset, width).trimmed(), value);
					QCOMPARE(line.mid(offset + (c < 2 ? 0 : width - value.size()), value.size()), value);
					offset += width;
					if (c + 1 < table.columnCount()) {
						QCOMPARE(line.mid(offset, 2), QString("  "));
						offset += 2;
					}
				}
				QCOMPARE(offset, line.size());
			}
		}
	}

	void initTestCase()
	{
		QVERIFY(m_settings.isValid());
		QSettings::setDefaultFormat(QSettings::IniFormat);
		QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_settings.path());
		registerQml();
	}

	void init()
	{
		QTest::failOnWarning(
		    QRegularExpression(".*(Unable to assign|ReferenceError|TypeError|Binding loop|Cannot assign|failed to "
		                       "load|Cannot open|recursive rearrange).*"));
	}

	void dataSourcesRejectMalformedFiles()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QVERIFY(QDir(directory.path()).mkpath("weights"));
		QFile file(directory.filePath("weights/vfcc17.rules"));
		QVERIFY(file.open(QIODevice::WriteOnly));
		file.write("not a valid weight table\n");
		file.close();
		const auto previous = getExecutableDir();
		const auto restore = qScopeGuard([&] { getExecutableDir() = previous; });
		getExecutableDir() = directory.path().toStdWString();
		DataSourcesModel sources;
		QVERIFY(!sources.hasWarning());
		sources.refresh();
		QTRY_VERIFY_WITH_TIMEOUT(!sources.checking(), 30000);
		QTRY_COMPARE_WITH_TIMEOUT(sources.entries().size(), 12, 30000);
		bool found = false, foundEfcc18 = false, foundFcct14 = false;
		for (const auto& value : sources.entries()) {
			const auto entry = value.toMap();
			const auto id = entry.value("id").toString();
			QVERIFY(id != "debias_2018");
			if (id == "debias_efcc18") {
				foundEfcc18 = true;
				QCOMPARE(entry.value("title").toString(), QString("Astrometric debias: EFCC18"));
			} else if (id == "debias_fcct14") {
				foundFcct14 = true;
				QCOMPARE(entry.value("title").toString(), QString("Astrometric debias: FCCT14"));
			}
			if (entry.value("id") != "weights")
				continue;
			found = true;
			QVERIFY(!entry.value("valid").toBool());
			QCOMPARE(entry.value("status").toString(), QString("Invalid"));
		}
		QVERIFY(found);
		QVERIFY(foundEfcc18);
		QVERIFY(foundFcct14);
		QVERIFY(sources.hasWarning());
	}

	void dataSourcesNoticeRemovedFile()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QVERIFY(QDir(directory.path()).mkpath("weights"));
		const QString weightPath = directory.filePath("weights/vfcc17.rules");
		QFile weights(weightPath);
		QVERIFY(weights.open(QIODevice::WriteOnly));
		weights.write("invalid weight data\n");
		weights.close();
		const auto previous = getExecutableDir();
		const auto restore = qScopeGuard([&] { getExecutableDir() = previous; });
		getExecutableDir() = directory.path().toStdWString();
		DataSourcesModel sources;
		sources.refresh();
		QTRY_VERIFY_WITH_TIMEOUT(sources.checked(), 30000);
		const auto weightEntry = [&] {
			for (const auto& value : sources.entries()) {
				const auto entry = value.toMap();
				if (entry.value("id") == "weights")
					return entry;
			}
			return QVariantMap{};
		};
		QVERIFY(weightEntry().value("onDisk").toBool());
		QVERIFY(QFile::remove(weightPath));
		QTRY_VERIFY_WITH_TIMEOUT(!weightEntry().value("onDisk").toBool(), 30000);
		QVERIFY(sources.restartRequired());
		QVERIFY(sources.unavailableCount() > 0);
	}

	void progressIndicatorStates()
	{
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import TrajectaUI
            Window {
                width: 1120; height: 112; visible: true; color: AppTheme.base
                property string theme: "graphite"
                Binding { target: AppTheme; property: "name"; value: theme }
                Item {
                    objectName: "progress-host"; anchors.fill: parent
                    CustomText {
                        x: 24; y: 24
                        text: "Preparing atmospheric entries and environment..."
                        font.pixelSize: 12; color: AppTheme.muted
                    }
                    CustomProgressBar {
                        objectName: "progress"; x: 24; y: 64
                        width: parent.width - 48; indeterminate: true
                    }
                }
            }
        )",
		                  QUrl());
		std::unique_ptr<QObject> root(component.create());
		QVERIFY2(root, qPrintable(component.errorString()));
		auto* window = qobject_cast<QQuickWindow*>(root.get());
		QVERIFY(window && QTest::qWaitForWindowExposed(window));
		auto* host = findItem(window->contentItem(), "progress-host");
		auto* progress = findItem(host, "progress");
		auto* indicator = findItem(progress, "progress-indicator");
		QVERIFY(progress && indicator);
		QTRY_VERIFY(progress->property("animating").toBool());
		const double start = indicator->x();
		const auto first = window->grabWindow();
		QTest::qWait(400);
		QVERIFY(std::abs(indicator->x() - start) > 5);
		QVERIFY(first != window->grabWindow());
		for (const auto& theme : {"graphite", "light", "dim"}) {
			window->setProperty("theme", theme);
			for (int width : {360, 1120, 2240}) {
				window->resize(width, 112);
				QTest::qWait(100);
				QVERIFY(indicator->x() >= 0);
				QVERIFY(indicator->x() + indicator->width() <= progress->width() + 0.01);
				QVERIFY(indicator->width() <= 180);
				QCOMPARE(
				    indicator->property("color"),
				    engine.singletonInstance<QObject*>(qmlTypeId("TrajectaUI", 1, 0, "AppTheme"))->property("accent"));
				const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
				if (!captures.isEmpty() && width == 1120) {
					QVERIFY(QDir().mkpath(captures));
					QVERIFY(window->grabWindow().save(captures + "/progress-" + theme + ".png"));
				}
			}
		}
		host->setVisible(false);
		QVERIFY(!progress->property("animating").toBool());
		QTest::qWait(50);
		const auto stopped = indicator->x();
		QTest::qWait(100);
		QCOMPARE(indicator->x(), stopped);
		host->setVisible(true);
		QTRY_VERIFY(progress->property("animating").toBool());
		window->hide();
		QVERIFY(!progress->property("animating").toBool());
		window->show();
		progress->setProperty("indeterminate", false);
		progress->setProperty("value", 0.42);
		QVERIFY(!progress->property("animating").toBool());
		QTRY_VERIFY(std::abs(indicator->width() - progress->width() * 0.42) < 0.01);
		QCOMPARE(indicator->x(), 0.0);
		progress->setProperty("value", 1.0);
		QTRY_VERIFY(std::abs(indicator->width() - progress->width()) < 0.01);
		progress->setProperty("value", 0.0);
		QTRY_VERIFY(indicator->width() < 0.01);
	}

	void inlineValueBaselines()
	{
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import TrajectaUI
            Window {
                width: 400; height: 180; visible: true
                InlineValue { objectName: "inline"; x: 20; y: 20; label: "q"; value: "0.98500507"; unit: "AU" }
                KeyValueList { objectName: "list"; x: 20; y: 60; width: 300; rows: [["Period", "431.4673 days"]] }
            }
        )",
		                  QUrl());
		std::unique_ptr<QObject> root(component.create());
		QVERIFY2(root, qPrintable(component.errorString()));
		auto* window = qobject_cast<QQuickWindow*>(root.get());
		auto* row = findItem(window->contentItem(), "inline");
		QVERIFY(row);
		for (int size : {11, 14, 18}) {
			row->setProperty("valueSize", size);
			QTest::qWait(50);
			const auto* label = findItem(row, "inline-label");
			const auto* value = findItem(row, "inline-value");
			const auto* unit = findItem(row, "inline-unit");
			QVERIFY(label && value && unit);
			const auto baseline = [](const QQuickItem* item) { return item->y() + item->baselineOffset(); };
			QVERIFY(std::abs(baseline(label) - baseline(value)) < 0.01);
			QVERIFY(std::abs(baseline(unit) - baseline(value)) < 0.01);
			QVERIFY(row->height() >= value->height());
		}
	}

	void riskDraftValidation()
	{
		SolutionViewModel solution;
		RiskMonitorViewModel monitor(solution);
		QVERIFY(!monitor.available());
		QVERIFY(!monitor.draft());
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(result, "risk-fixture"));
		QVERIFY(monitor.available());
		QCOMPARE(monitor.covarianceSize(), QStringLiteral("6 x 6"));
		QCOMPARE(monitor.observationCount(), 2);
		QVERIFY(monitor.sourceError().isEmpty());
		monitor.setCenterUtc("2023-02-27 00:00:00");
		QVERIFY2(monitor.inputReady(), qPrintable(monitor.formError()));
		const auto original = monitor.draft().value();
		QCOMPARE(original.endMjdUtc - original.startMjdUtc, 2.0);
		QCOMPARE(original.startMjdUtc, 60001.0);
		monitor.setSeed("18446744073709551615");
		QCOMPARE(monitor.draft()->seed, std::numeric_limits<std::uint64_t>::max());
		for (const auto* seed : {"18446744073709551616", "-1", "1e3", "1.0", "", " 2"}) {
			monitor.setSeed(QString::fromLatin1(seed));
			QVERIFY(!monitor.draft());
		}
		monitor.setSeed("0");
		QCOMPARE(monitor.draft()->seed, std::uint64_t{0});
		for (const auto* count : {"0", "-1", "1.2", "1e3", "2147483648", ""}) {
			monitor.setSamples(QString::fromLatin1(count));
			QVERIFY(!monitor.draft());
		}
		monitor.setSamples("2000");
		monitor.setMethod(1);
		QCOMPARE(monitor.samples(), QStringLiteral("1000"));
		monitor.setSamples("50");
		monitor.setMethod(0);
		QCOMPARE(monitor.samples(), QStringLiteral("2000"));
		monitor.setMethod(99);
		QCOMPARE(monitor.method(), 0);
		monitor.setTarget(1);
		QCOMPARE(monitor.boundary(), 0);
		QCOMPARE(monitor.altitudeKm(), 0.0);
		monitor.setTarget(0);
		QCOMPARE(monitor.boundary(), 1);
		QCOMPARE(monitor.altitudeKm(), 100.0);
		monitor.setBoundary(0);
		QVERIFY(monitor.draft());
		monitor.setCenterUtc("2023-02-25 00:00:00");
		QVERIFY(!monitor.draft());
		monitor.setCenterUtc("2023-02-27 00:00:00");
		monitor.setHalfWindowHours(0);
		QVERIFY(!monitor.draft());
		monitor.setHalfWindowHours(std::numeric_limits<double>::max());
		QVERIFY(!monitor.draft());
		monitor.setHalfWindowHours(std::numeric_limits<double>::quiet_NaN());
		QVERIFY(!monitor.draft());
		monitor.setHalfWindowHours(24);
		monitor.setCenterUtc("2023-02-30 00:00:00");
		QVERIFY(!monitor.draft());
		monitor.setCenterUtc("2023-02-27 24:00:00");
		QVERIFY(!monitor.draft());
		monitor.setCenterUtc("2023-02-27 00:00:00");
		QVERIFY(solution.publish(result, "second-source"));
		QVERIFY(monitor.draft()->sourceRevision > original.sourceRevision);
		QCOMPARE(monitor.centerUtc(), QStringLiteral("2023-02-27 00:00:00"));
		QCOMPARE(monitor.samples(), QStringLiteral("2000"));
		result.fit.uncertainty.fitParameterCovariance(0, 0) = -1.0;
		QVERIFY(solution.publish(result, "bad-covariance"));
		QVERIFY(!monitor.draft());
		monitor.setMethod(1);
		QVERIFY(monitor.draft());
		QCOMPARE(monitor.draft()->samples, 50);
		result.fit.finalNonGrav.setModel(NgModel::Yarkovsky);
		result.fit.finalNonGrav.setEnabled(true);
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(6, 6) * 1e-8;
		QVERIFY(solution.publish(result, "non-grav-state-only"));
		monitor.setMethod(0);
		QVERIFY(!monitor.draft());
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(7, 7) * 1e-8;
		QVERIFY(solution.publish(result, "full-non-grav"));
		QVERIFY(monitor.draft());
		QCOMPARE(monitor.covarianceSize(), QStringLiteral("7 x 7"));
		result.fit.uncertainty.fitParameterCovariance(0, 1) = 2e-8;
		result.fit.uncertainty.fitParameterCovariance(1, 0) = 2e-8;
		QVERIFY(solution.publish(result, "indefinite-non-grav"));
		QVERIFY(!monitor.draft());
		monitor.setMethod(1);
		result.fit.converged = false;
		QVERIFY(solution.publish(result, "unconverged"));
		QVERIFY(!monitor.draft());
		result.fit.fitEpochOrbit.reset();
		QVERIFY(solution.publish(result, "no-fit-epoch"));
		QVERIFY(!monitor.available());
		QVERIFY(!monitor.draft());
	}

	void riskBoundaryPresets()
	{
		SolutionViewModel solution;
		RiskMonitorViewModel monitor(solution);
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(result, "risk-boundary"));
		monitor.setCenterUtc("2023-02-27 00:00:00");
		const auto altitude = monitor.metaObject()->property(monitor.metaObject()->indexOfProperty("altitudeKm"));
		QVERIFY(altitude.isValid());
		QVERIFY(!altitude.isWritable());
		for (int method : {0, 1}) {
			monitor.setMethod(method);
			for (int interval : {0, 1}) {
				monitor.setIntervalMode(interval);
				monitor.setTarget(0);
				monitor.setBoundary(1);
				QVERIFY(monitor.draft());
				QCOMPARE(monitor.draft()->altitudeKm, 100.0);
				QVERIFY(!monitor.setProperty("altitudeKm", 20.0));
				QCOMPARE(monitor.draft()->altitudeKm, 100.0);
				for (int target : {1, 2}) {
					monitor.setTarget(target);
					monitor.setBoundary(1);
					QCOMPARE(monitor.boundary(), 0);
					QCOMPARE(monitor.draft()->altitudeKm, 0.0);
				}
				monitor.setTarget(0);
				QCOMPARE(monitor.boundary(), 1);
				QCOMPARE(monitor.draft()->altitudeKm, 100.0);
				monitor.setBoundary(0);
				QCOMPARE(monitor.draft()->altitudeKm, 0.0);
				monitor.setBoundary(99);
				QCOMPARE(monitor.draft()->altitudeKm, 0.0);
			}
		}
	}

	void atmosphereDraftValidation()
	{
		SolutionViewModel solution;
		RiskMonitorViewModel risk(solution);
		AtmosphereViewModel atmosphere(solution, risk);
		QVERIFY(!atmosphere.sourceError().isEmpty());
		QVERIFY(!atmosphere.formError().isEmpty());
		QCOMPARE(atmosphere.preset(), 0);
		atmosphere.setValue("magnitude", "25.3");
		QVERIFY2(atmosphere.formError().isEmpty(), qPrintable(atmosphere.formError()));
		const auto original = atmosphere.parameterOptions();
		QCOMPARE(original.physics.ramPressure, strewn::RamPressure::RhoV2);
		QCOMPARE(original.physics.representativesPerBin, std::size_t{3});
		QCOMPARE(original.monteCarlo.parameters.weakProbability, 0.0);
		for (int preset = 0; preset < 3; ++preset) {
			atmosphere.selectPreset(preset);
			QCOMPARE(atmosphere.preset(), preset);
			QCOMPARE(atmosphere.values()["magnitude"].toString(), QStringLiteral("25.3"));
			QVERIFY(atmosphere.formError().isEmpty());
		}
		atmosphere.selectPreset(0);
		atmosphere.setValue("strengthMin", "6");
		QCOMPARE(atmosphere.preset(), 3);
		QVERIFY(!atmosphere.formError().isEmpty());
		atmosphere.selectPreset(0);
		atmosphere.setValue("seed", "18446744073709551615");
		QCOMPARE(atmosphere.parameterOptions().monteCarlo.seed, std::numeric_limits<std::uint64_t>::max());
		for (const auto* seed : {"18446744073709551616", "-1", "1e3", "1.0", "", " 2"}) {
			atmosphere.setValue("seed", QString::fromLatin1(seed));
			QVERIFY(!atmosphere.formError().isEmpty());
		}
		atmosphere.setValue("seed", "42");
		atmosphere.setValue("massMin", "200");
		QVERIFY(atmosphere.formError().isEmpty());
		atmosphere.setSizeMode(1);
		QVERIFY(!atmosphere.formError().isEmpty());
		atmosphere.setValue("massMax", "300");
		atmosphere.setValue("magnitude", "");
		QVERIFY(atmosphere.formError().isEmpty());
		QCOMPARE(std::get<strewn::MassSize>(atmosphere.parameterOptions().monteCarlo.parameters.size).kg.min, 200.0);
		atmosphere.selectPreset(2);
		QCOMPARE(atmosphere.sizeMode(), 1);
		QCOMPARE(atmosphere.values()["massMin"].toString(), QStringLiteral("200"));
		atmosphere.setSizeMode(2);
		QVERIFY(atmosphere.formError().isEmpty());
		atmosphere.setSourceMode(1);
		atmosphere.setValue("epoch", "2008-10-07 02:45:00");
		QVERIFY(atmosphere.sourceError().isEmpty());
		for (const auto* epoch : {"2008-02-30 00:00:00", "2008-10-07 24:00:00", ""}) {
			atmosphere.setValue("epoch", QString::fromLatin1(epoch));
			QVERIFY(!atmosphere.sourceError().isEmpty());
		}
		atmosphere.setValue("epoch", "2008-10-07 02:45:00");
		atmosphere.setValue("slope", "1");
		QVERIFY(!atmosphere.sourceError().isEmpty());
		atmosphere.setValue("slope", "-90");
		QVERIFY(atmosphere.sourceError().isEmpty());
		atmosphere.setValue("threshold", "");
		QVERIFY(!atmosphere.formError().isEmpty());
		atmosphere.setConvergence(false);
		QVERIFY(atmosphere.formError().isEmpty());
		QVERIFY(!atmosphere.parameterOptions().monteCarlo.convergence);
		atmosphere.setOffline(true);
		QVERIFY(std::get<strewn::AutomaticEnvironment>(atmosphere.parameterOptions().environment).offline);
		atmosphere.setEnvironmentMode(2);
		QVERIFY(atmosphere.formError().isEmpty());
		atmosphere.setValue("workers", "0");
		QVERIFY(!atmosphere.formError().isEmpty());
		const auto values = atmosphere.values();
		atmosphere.setValue("height", "20");
		QCOMPARE(atmosphere.values(), values);
	}

	void atmosphereWorkspaceLayout()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->setProperty("currentPage", 4);
		QTest::qWait(80);
		auto* page = findItem(window->contentItem(), "risk-monitor-page");
		QVERIFY(page);
		auto* tab = findItem(page, "risk-method-2");
		QVERIFY(tab);
		QVERIFY(QMetaObject::invokeMethod(tab, "clicked"));
		QTest::qWait(80);
		auto* form = findItem(page, "atmosphere-form");
		QVERIFY(form && form->isVisible());
		QVERIFY(!findItem(page, "risk-center")->isVisible());
		QVERIFY(!findItem(page, "risk-run")->isEnabled());
		QVERIFY(!findItem(page, "atmosphere-height"));
		session.atmosphere()->setValue("magnitude", "25.3");
		session.atmosphere()->setValue("seed", "18446744073709551615");
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		auto* scroll = findItem(page, "risk-scroll");
		QVERIFY(scroll);
		auto* scrollBar = findItem(page, "risk-scrollbar");
		QVERIFY(scrollBar);
		auto* flickable = qvariant_cast<QObject*>(scroll->property("contentItem"));
		QVERIFY(flickable);
		for (const auto* theme : {"graphite", "light"}) {
			session.setTheme(QString::fromLatin1(theme));
			for (const QSize size : {QSize(1440, 900), QSize(1920, 1080), QSize(1024, 720)}) {
				window->resize(size);
				for (int source : {0, 1}) {
					session.atmosphere()->setSourceMode(source);
					QTest::qWait(100);
					QVERIFY(scrollBar->isVisible());
					QCOMPARE(scrollBar->parentItem(), scroll);
					QCOMPARE(scrollBar->height(), scroll->property("availableHeight").toDouble());
					QCOMPARE(scrollBar->x() + scrollBar->width(), scroll->width());
					QCOMPARE(scrollBar->property("orientation").toInt(), int(Qt::Vertical));
					QVERIFY(scrollBar->property("size").toDouble() < 1.0);
					flickable->setProperty("contentY", 0);
					for (const auto* key : {"densityMin", "densityMax", "strengthMin", "albedoMin", "albedoMode",
					                        "albedoMax", "seed", "runs"}) {
						auto* field = findItem(page, "atmosphere-" + QString::fromLatin1(key));
						QVERIFY(field && field->isVisible());
						const auto rect = field->mapRectToItem(form, field->boundingRect());
						QVERIFY2(rect.left() >= 0 && rect.right() <= form->width() + 1, key);
						QVERIFY2(field->width() >= 65, key);
					}
					QCOMPARE(findItem(page, "atmosphere-seed")->property("text").toString(),
					         QStringLiteral("18446744073709551615"));
					for (const auto* name : {"atmosphere-source-mode", "atmosphere-preset", "atmosphere-size-mode",
					                         "atmosphere-environment-mode", "atmosphere-map-method"}) {
						auto* field = findItem(page, QString::fromLatin1(name));
						QVERIFY(field && field->isVisible());
						QVERIFY2(field->width() >= 200, name);
						const auto rect = field->mapRectToItem(form, field->boundingRect());
						QVERIFY2(rect.left() >= 0 && rect.right() <= form->width() + 1, name);
					}
					const auto prefix = captures + "/atmosphere-" + theme + "-" + QString::number(size.width()) + "-" +
					                    QString::number(source);
					QTest::qWait(30);
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(prefix + "-top.png"));
					flickable->setProperty("contentY", 700);
					QTest::qWait(30);
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(prefix + "-material.png"));
					const double bottom =
					    flickable->property("contentHeight").toDouble() - flickable->property("height").toDouble();
					QVERIFY(bottom > 0);
					flickable->setProperty("contentY", std::max(0.0, bottom));
					QTest::qWait(30);
					QVERIFY(flickable->property("contentY").toDouble() > 0);
					QVERIFY(scrollBar->property("position").toDouble() > 0);
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(prefix + "-bottom.png"));
				}
			}
		}
		auto* density = findItem(page, "atmosphere-densityMin");
		QVERIFY(density->setProperty("text", "1000"));
		QVERIFY(QMetaObject::invokeMethod(density, "textEdited"));
		QCOMPARE(session.atmosphere()->values()["densityMin"].toString(), QStringLiteral("1000"));
		QCOMPARE(session.atmosphere()->preset(), 3);
		auto* preset = findItem(page, "atmosphere-preset");
		QVERIFY(preset->setProperty("currentIndex", 2));
		QVERIFY(QMetaObject::invokeMethod(preset, "activated", Q_ARG(int, 2)));
		QCOMPARE(session.atmosphere()->preset(), 2);
		QCOMPARE(density->property("text").toString(), QStringLiteral("7000"));
		QVERIFY(QMetaObject::invokeMethod(findItem(page, "risk-method-0"), "clicked"));
		QVERIFY(!findItem(page, "atmosphere-form")->isVisible());
		QVERIFY(QMetaObject::invokeMethod(tab, "clicked"));
		QCOMPARE(session.atmosphere()->values()["seed"].toString(), QStringLiteral("18446744073709551615"));
		QVERIFY(!findItem(page, "risk-run")->isEnabled());
	}

	void closeApproachesScrollbarGeometry()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->resize(1024, 720);
		window->setProperty("currentPage", 5);
		QTest::qWait(100);
		auto* page = findItem(window->contentItem(), "close-approaches-page");
		QVERIFY(page && page->isVisible());
		auto* scroll = findItem(page, "approaches-scroll");
		auto* scrollBar = findItem(page, "approaches-scrollbar");
		QVERIFY(scroll && scrollBar);
		QVERIFY(scrollBar->isVisible());
		QCOMPARE(scrollBar->parentItem(), scroll);
		QCOMPARE(scrollBar->height(), scroll->property("availableHeight").toDouble());
		QCOMPARE(scrollBar->x() + scrollBar->width(), scroll->width());
		QCOMPARE(scrollBar->property("orientation").toInt(), int(Qt::Vertical));
		auto* flickable = qvariant_cast<QObject*>(scroll->property("contentItem"));
		QVERIFY(flickable);
		const double bottom =
		    flickable->property("contentHeight").toDouble() - flickable->property("height").toDouble();
		QVERIFY(bottom > 0);
		flickable->setProperty("contentY", bottom);
		QTest::qWait(30);
		QVERIFY(scrollBar->property("position").toDouble() > 0);
	}

	void atmosphereArchiveValidation()
	{
		QTemporaryDir temp;
		const auto path = temp.filePath("result.json");
		auto document = atmosphericArchiveFixture();
		atmosphere::saveJson(path, document);
		const auto archive = atmosphere::readArchive(path);
		QCOMPARE(archive.events.size(), std::size_t(1));
		QCOMPARE(archive.events[0].entries.size(), std::size_t(1));
		const auto& e = archive.events[0].entries[0];
		QCOMPARE(e.id, std::uint64_t(0));
		QCOMPARE(e.weight, 1.0);
		const auto expected = strewn::enuToEcef({0, 500, -1000}, 0.5, 0.3);
		for (int i = 0; i < 3; ++i)
			QVERIFY(std::abs(e.velocityEcefMps[i] - expected[i]) < 1e-6);
		QVERIFY(std::abs(strewn::geodetic(e.positionEcefM).heightWgs84M - 100000) < 1e-6);
		QVERIFY(std::abs(orbital_core::time::toMjdTt(e.epochMjdUtc, orbital_core::time::Scale::UTC) - 60000.25) <
		        1e-10);
		atmosphere::verifyArchive(archive);
		for (const auto& change : std::vector<std::pair<QString, QJsonValue>>{{"schema", "bad"},
		                                                                      {"time_scale", "MJD UTC"},
		                                                                      {"target", 8},
		                                                                      {"altitude_km", 20},
		                                                                      {"cancelled", true},
		                                                                      {"completed", 1},
		                                                                      {"failed_fit", 1},
		                                                                      {"entry_clones", 2},
		                                                                      {"body_fixed_frame", "unknown"},
		                                                                      {"body_fixed_model", "unknown"}}) {
			auto bad = document;
			bad[change.first] = change.second;
			atmosphere::saveJson(path, bad);
			QVERIFY_EXCEPTION_THROWN(atmosphere::readArchive(path), std::exception);
		}
		auto rows = document["clones"].toArray();
		rows[1] = rows[0];
		auto bad = document;
		bad["clones"] = rows;
		atmosphere::saveJson(path, bad);
		QVERIFY_EXCEPTION_THROWN(atmosphere::readArchive(path), std::exception);
		QVERIFY_EXCEPTION_THROWN(atmosphere::verifyArchive(archive), std::exception);
		document["method"] = "observations";
		atmosphere::saveJson(path, document);
		QCOMPARE(atmosphere::readArchive(path).events[0].entries.size(), std::size_t(1));
		const auto config = strewn::atmosphericConvergenceGrid(archive.events[0].entries);
		QVERIFY(std::abs(config.centerDeg[0] - 0.5 * 180 / Constants::PI) < 1);
		QVERIFY(config.halfSizeKm[0] > 0);
	}

	void atmosphereArchiveDoesNotFollowAnotherObject()
	{
		QTemporaryDir temp;
		const auto path = temp.filePath("result.json");
		atmosphere::saveJson(path, atmosphericArchiveFixture());
		SolutionViewModel solution;
		RiskMonitorViewModel risk(solution);
		AtmosphereViewModel atmosphere(solution, risk);
		auto first = fixture();
		first.observations.objectName = "First";
		QVERIFY(solution.publish(first, "first-object"));
		atmosphere.selectArchive(QUrl::fromLocalFile(path));
		QTRY_VERIFY_WITH_TIMEOUT(!atmosphere.busy(), 10000);
		QVERIFY2(atmosphere.sourceError().isEmpty(), qPrintable(atmosphere.sourceError()));
		QVERIFY(solution.publish(first, "first-object-refit"));
		QCOMPARE(atmosphere.archivePath(), path);
		auto second = fixture();
		second.observations.objectName = "Second";
		QVERIFY(solution.publish(second, "second-object"));
		QVERIFY(atmosphere.archivePath().isEmpty());
		QCOMPARE(atmosphere.selectedEvent(), -1);
		QVERIFY(!atmosphere.sourceError().isEmpty());

		atmosphere.selectArchive(QUrl::fromLocalFile(path));
		QVERIFY(solution.publish(first, "first-object-again"));
		QTRY_VERIFY_WITH_TIMEOUT(!atmosphere.busy(), 10000);
		QVERIFY(atmosphere.archivePath().isEmpty());
		QCOMPARE(atmosphere.selectedEvent(), -1);
	}

	void atmosphereRunLifecycle()
	{
		QTemporaryDir temp;
		const auto path = temp.filePath("result.json");
		atmosphere::saveJson(path, atmosphericArchiveFixture());
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto* model = session.atmosphere();
		model->selectArchive(QUrl::fromLocalFile(path));
		QVERIFY(session.busy());
		QTRY_VERIFY_WITH_TIMEOUT(!model->busy(), 10000);
		QVERIFY2(model->sourceError().isEmpty(), qPrintable(model->sourceError()));
		QCOMPARE(model->selectedEvent(), 0);
		QVERIFY(model->values()["magnitude"].toString().isEmpty());
		model->setSizeMode(1);
		model->setValue("massMin", "0.01");
		model->setValue("massMax", "0.01");
		model->setValue("runs", "2");
		model->setValue("workers", "2");
		model->setEnvironmentMode(2);
		model->setOffline(true);
		model->setConvergence(false);
		QVERIFY2(model->canRun(), qPrintable(model->sourceError() + model->formError()));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->setProperty("currentPage", 4);
		QTest::qWait(80);
		auto* page = findItem(window->contentItem(), "risk-monitor-page");
		QVERIFY(page);
		page->setProperty("showingAtmosphere", true);
		auto* run = findItem(page, "risk-run");
		QVERIFY(run && run->isEnabled());
		QVERIFY(QMetaObject::invokeMethod(run, "clicked"));
		QCOMPARE(page->property("currentView").toInt(), 1);
		QVERIFY(session.busy());
		QVERIFY(!session.risk()->canRun());
		QVERIFY(!session.canRun());
		const auto values = model->values();
		model->setValue("runs", "200");
		QCOMPARE(model->values(), values);
		QTRY_VERIFY_WITH_TIMEOUT(!model->busy(), 60000);
		QCOMPARE(model->runs().size(), 1);
		QVERIFY2(model->error().isEmpty(), qPrintable(model->error()));
		QVERIFY2(model->reportAvailable(), qPrintable(model->summary()));
		QTRY_VERIFY(findItem(page, "atmosphere-run-view")->isVisible());
		QVERIFY(!findItem(page, "risk-runs")->isVisible());
		QVERIFY(findItem(page, "atmosphere-report")->isEnabled());
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty()) {
			QVERIFY(QDir().mkpath(captures));
			QTest::qWait(80);
			QVERIFY(window->grabWindow().save(captures + "/atmosphere-result.png"));
		}
		const auto successful = model->resultPath();
		QVERIFY(QFileInfo(successful + "/launch.json").isFile());
		QVERIFY(QFileInfo(successful + "/source.json").isFile());
		QVERIFY(QFileInfo(successful + "/results/summary.json").isFile());
		const auto replay = strewn::loadJob(std::filesystem::path((successful + "/launch.json").toStdWString()));
		QCOMPARE(replay.monteCarlo.runs, std::size_t(2));
		QCOMPARE(std::get<strewn::IdealizedEnvironment>(replay.environment).groundHeightWgs84M, 0.0);
		model->run();
		session.cancel();
		QVERIFY(model->cancelling());
		QTRY_VERIFY_WITH_TIMEOUT(!model->busy(), 30000);
		QCOMPARE(model->runs().size(), 2);
		QVERIFY(model->summary().contains("Cancelled"));
		model->setSelectedRun(0);
		QCOMPARE(model->resultPath(), successful);
		QVERIFY(model->reportAvailable());
		auto changed = atmosphericArchiveFixture();
		changed["seed"] = "changed";
		atmosphere::saveJson(path, changed);
		model->run();
		QTRY_VERIFY_WITH_TIMEOUT(!model->busy(), 10000);
		QVERIFY(model->error().contains("changed after selection"));
		QVERIFY(model->resultPath().isEmpty());
		model->setSelectedRun(0);
		QVERIFY(model->reportAvailable());
	}

	void atmosphereSavedArchiveRun()
	{
		const auto source = qEnvironmentVariable("ORBITAL_ATMOSPHERE_ARCHIVE");
		if (source.isEmpty())
			QSKIP("Saved orbital archive integration run is opt-in");
		const auto archive = atmosphere::readArchive(source);
		QVERIFY(!archive.events.empty());
		QVERIFY(archive.magnitude.has_value());
		SolutionViewModel solution;
		RiskMonitorViewModel risk(solution);
		AtmosphereViewModel model(solution, risk);
		model.selectArchive(QUrl::fromLocalFile(source));
		QTRY_VERIFY_WITH_TIMEOUT(!model.busy(), 30000);
		model.setSelectedEvent(0);
		QCOMPARE(model.values()["magnitude"].toString().toDouble(), *archive.magnitude);
		model.setValue("runs", "2");
		model.setValue("workers", "2");
		model.setValue("timeout", "60");
		model.setEnvironmentMode(qEnvironmentVariableIsSet("ORBITAL_ATMOSPHERE_AUTO") ? 0 : 2);
		model.setOffline(model.environmentMode() != 0);
		model.setConvergence(false);
		QVERIFY2(model.canRun(), qPrintable(model.sourceError() + model.formError()));
		model.run();
		QTRY_VERIFY_WITH_TIMEOUT(!model.busy(), 300000);
		qInfo().noquote() << model.resultPath() << model.summary() << model.error();
		QVERIFY2(model.error().isEmpty(), qPrintable(model.error()));
		QVERIFY(model.reportAvailable());
		atmosphere::verifyArchive(archive);
		QFile saved(model.resultPath() + "/entries.json");
		QVERIFY(saved.open(QIODevice::ReadOnly));
		const auto rows = QJsonDocument::fromJson(saved.readAll()).object()["entries"].toArray();
		QCOMPARE(rows.size(), qsizetype(archive.events[0].entries.size()));
	}

	void riskWorkspaceLayout()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->setProperty("currentPage", 4);
		QTest::qWait(100);
		auto* page = findItem(window->contentItem(), "risk-monitor-page");
		QVERIFY(page);
		auto* run = findItem(page, "risk-run");
		QVERIFY(page && page->isVisible() && run && !run->isEnabled());
		QVERIFY(!session.risk()->inputReady());
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty()) {
			QVERIFY(QDir().mkpath(captures));
			QVERIFY(window->grabWindow().save(captures + "/risk-empty.png"));
		}
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		result.observations.objectName = "K20V00V";
		QVERIFY(session.solution()->publish(result, "risk-layout"));
		auto* center = findItem(page, "risk-center");
		QVERIFY(center && !center->property("text").toString().isEmpty());
		QCOMPARE(center->property("text").toString(), session.risk()->centerUtc());
		QVERIFY(!findItem(page, "risk-window-start"));
		QVERIFY(!findItem(page, "risk-altitude"));
		QVERIFY(!findItem(page, "risk-range-start")->isVisible());
		QVERIFY(session.risk()->inputReady());
		session.risk()->setCenterUtc("2044-10-10 12:00:00");
		session.risk()->setSeed("18446744073709551615");
		for (const auto* theme : {"graphite", "light"}) {
			session.setTheme(QString::fromLatin1(theme));
			for (const QSize size : {QSize(1440, 900), QSize(1440, 1080), QSize(1024, 720)}) {
				window->resize(size);
				for (int method = 0; method < 2; ++method) {
					session.risk()->setMethod(method);
					session.risk()->setIntervalMode(0);
					QTest::qWait(120);
					QVERIFY(session.risk()->inputReady());
					QVERIFY(run->isEnabled());
					auto* boundary = findItem(page, "risk-boundary");
					QVERIFY(boundary);
					QCOMPARE(boundary->property("count").toInt(), 2);
					QCOMPARE(boundary->property("currentText").toString(),
					         QStringLiteral("Atmospheric entry \u00b7 100 km"));
					QVERIFY(boundary->setProperty("currentIndex", 0));
					QVERIFY(QMetaObject::invokeMethod(boundary, "activated", Q_ARG(int, 0)));
					QCOMPARE(session.risk()->draft()->altitudeKm, 0.0);
					QVERIFY(boundary->setProperty("currentIndex", 1));
					QVERIFY(QMetaObject::invokeMethod(boundary, "activated", Q_ARG(int, 1)));
					QCOMPARE(session.risk()->draft()->altitudeKm, 100.0);
					for (const auto* name : {"risk-center", "risk-half-window", "risk-target", "risk-boundary",
					                         "risk-samples", "risk-seed"}) {
						auto* control = findItem(page, QString::fromLatin1(name));
						QVERIFY(control && control->isVisible());
						const auto rect = control->mapRectToItem(page, control->boundingRect());
						QVERIFY2(rect.left() >= 0 && rect.right() <= page->width(), name);
						QVERIFY2(control->width() >= 200, name);
					}
					QCOMPARE(findItem(page, "risk-seed")->property("text").toString(),
					         QStringLiteral("18446744073709551615"));
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(
						    captures + "/risk-" + theme + "-" + QString::number(size.width()) + "x" +
						    QString::number(size.height()) + "-" + QString::number(method) + ".png"));
					auto* mode = findItem(page, "risk-date-source");
					QVERIFY(mode);
					QVERIFY(mode->setProperty("currentIndex", 1));
					QVERIFY(QMetaObject::invokeMethod(mode, "activated", Q_ARG(int, 1)));
					QCOMPARE(session.risk()->intervalMode(), 1);
					QTest::qWait(120);
					QVERIFY2(session.risk()->inputReady(), qPrintable(session.risk()->formError()));
					QVERIFY(!findItem(page, "risk-center")->isVisible());
					QVERIFY(!findItem(page, "risk-half-window")->isVisible());
					for (const auto* name : {"risk-range-start", "risk-range-end", "risk-range-years"}) {
						auto* control = findItem(page, QString::fromLatin1(name));
						QVERIFY(control && control->isVisible());
						const auto rect = control->mapRectToItem(page, control->boundingRect());
						QVERIFY2(rect.left() >= 0 && rect.right() <= page->width(), name);
						QVERIFY2(control->width() >= 200, name);
					}
					QVERIFY(QMetaObject::invokeMethod(findItem(page, "risk-range-years"), "valueEdited",
					                                  Q_ARG(double, 10.0)));
					QCOMPARE(session.risk()->rangeYears(), 10.0);
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(
						    captures + "/risk-range-" + theme + "-" + QString::number(size.width()) + "x" +
						    QString::number(size.height()) + "-" + QString::number(method) + ".png"));
				}
			}
		}
		const auto sourceRevision = session.risk()->draft()->sourceRevision;
		page->setProperty("currentView", 1);
		QTest::qWait(50);
		QVERIFY(findItem(page, "risk-runs")->isVisible());
		session.risk()->events()->setEvents({{.startMjdTdb = 64001.1,
		                                      .endMjdTdb = 64001.2,
		                                      .entryClones = 12,
		                                      .probability = .0012,
		                                      .probabilityLow = .00068,
		                                      .probabilityHigh = .0021},
		                                     {.startMjdTdb = 65002.1,
		                                      .endMjdTdb = 65002.2,
		                                      .entryClones = 1,
		                                      .probability = .0001,
		                                      .probabilityHigh = .00057}});
		for (const auto* theme : {"graphite", "light"}) {
			session.setTheme(QString::fromLatin1(theme));
			for (const auto size : {QSize(1440, 900), QSize(1024, 720)}) {
				window->resize(size);
				QTest::qWait(100);
				QVERIFY(findItem(page, "risk-event-table")->isVisible());
				if (!captures.isEmpty())
					QVERIFY(window->grabWindow().save(captures + "/risk-events-" + theme + "-" +
					                                  QString::number(size.width()) + ".png"));
			}
		}
		QCOMPARE(session.risk()->draft()->sourceRevision, sourceRevision);
		page->setProperty("currentView", 0);
		QTest::qWait(50);
		QVERIFY(!findItem(page, "risk-runs")->isVisible());
		QVERIFY(QMetaObject::invokeMethod(findItem(page, "risk-method-0"), "clicked"));
		QCOMPARE(session.risk()->method(), 0);
		QVERIFY(QMetaObject::invokeMethod(findItem(page, "risk-method-1"), "clicked"));
		QCOMPARE(session.risk()->method(), 1);
		session.solution()->observations()->setEditingEnabled(true);
		session.solution()->observations()->setIncluded(1, true);
		QVERIFY(session.risk()->draft()->sourceRevision == sourceRevision);
		QCOMPARE(session.risk()->observationCount(), 2);
		QCOMPARE(session.solution()->snapshot()->fit.keptObsIndices.size(), std::size_t{2});
		QCOMPARE(result.fit.keptObsIndices.size(), std::size_t{2});
	}

	void riskJobLifecycle()
	{
		class Ephemeris final : public IEphemerisSource {
		   public:
			EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, {1, 0, 0}, {}}; }

			std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

			bool isAvailable() const override { return true; }

			std::string name() const override { return "risk-ui-test"; }
		} ephemeris;

		auto* old = getDefaultEphemerisSource();
		setDefaultEphemerisSource(&ephemeris);
		const auto restore = qScopeGuard([old] { setDefaultEphemerisSource(old); });
		QTemporaryDir output;
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		result.effectiveConfig.nbody.planetMask = 0;
		result.effectiveConfig.nbody.useSb441 = false;
		const auto residuals = result.fit.allResiduals.data();
		result.fit.allResiduals.clear();
		result.fit.residuals.clear();
		for (std::size_t i = 0; i < residuals.size(); ++i) {
			auto residual = residuals[i];
			residual.astrometricWeight = ResidualWeight{};
			result.fit.allResiduals.add(residual);
			if (i != 1)
				result.fit.residuals.add(residual);
		}
		QVERIFY(session.solution()->publish(result, output.path()));
		auto* monitor = session.risk();
		QAbstractItemModelTester tester(monitor->events(), QAbstractItemModelTester::FailureReportingMode::QtTest);
		monitor->setIntervalMode(1);
		monitor->setRangeYears(.001);
		monitor->setSamples("4");
		QVERIFY(monitor->canRun());
		monitor->run();
		QVERIFY(monitor->busy());
		QVERIFY(session.busy());
		QVERIFY(!session.canRun());
		QVERIFY(!session.ephemerides()->canGenerate());
		monitor->run();
		QTRY_VERIFY_WITH_TIMEOUT(!monitor->busy(), 20000);
		QVERIFY2(monitor->error().isEmpty(), qPrintable(monitor->error()));
		QCOMPARE(monitor->runs().size(), 1);
		QVERIFY(QFileInfo::exists(monitor->archivePath() + "/result.json"));
		QVERIFY(QFileInfo::exists(monitor->archivePath() + "/report.html"));
		const auto firstArchive = monitor->archivePath();
		monitor->setSamples("10000");
		monitor->run();
		monitor->cancel();
		QVERIFY(monitor->cancelling());
		QTRY_VERIFY_WITH_TIMEOUT(!monitor->busy(), 20000);
		QCOMPARE(monitor->runs().size(), 2);
		QVERIFY(monitor->warning().contains("Cancelled"));
		monitor->setSelectedRun(0);
		QCOMPARE(monitor->archivePath(), firstArchive);
		QVERIFY(session.solution()->publish(result, output.path()));
		QVERIFY(monitor->runs().isEmpty());
		QCOMPARE(monitor->events()->count(), 0);
	}

	void movValidationAndSourceChanges()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto* mov = session.mov();
		QVERIFY(!mov->canRun());
		auto source = fixture();
		source.effectiveConfig.nbody.enabled = true;
		source.effectiveConfig.useWeights = true;
		QVERIFY(session.solution()->publish(source, m_settings.path()));
		QVERIFY(mov->canRun());
		QVERIFY(!mov->seedUtc().isEmpty());
		mov->run({{"startUtc", "invalid"}});
		QVERIFY(!mov->error().isEmpty());
		QVERIFY(!mov->busy());
		QVERIFY(mov->runs().isEmpty());
		MovViewModel blocked(*session.solution(), [] { return OptimizationSettings{}; });
		blocked.setBlocked(true);
		QVERIFY(!blocked.canRun());
		blocked.setBlocked(false);
		QVERIFY(blocked.canRun());
		QVERIFY(session.solution()->publish(source, m_settings.path()));
		QVERIFY(mov->error().isEmpty());
	}

	void movFailureReturnsToUnchangedForm()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto source = fixture();
		source.effectiveConfig.nbody.enabled = true;
		source.effectiveConfig.useWeights = true;
		QVERIFY(session.solution()->publish(source, m_settings.path()));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->setProperty("currentPage", 4);
		auto* page = findItem(window->contentItem(), "risk-monitor-page");
		auto* form = findItem(window->contentItem(), "mov-form");
		QVERIFY(page);
		QVERIFY(form);
		page->setProperty("showingMov", true);
		const QString start = "2023-03-01 01:23:45", end = "2023-03-02 01:23:45";
		form->setProperty("startUtc", start);
		form->setProperty("endUtc", end);
		form->setProperty("sampling", 2);
		form->setProperty("observationSelection", 1);
		form->setProperty("rangeNodes", 16);
		form->setProperty("rateNodes", 24);
		auto* previous = getDefaultEphemerisSource();
		setDefaultEphemerisSource(nullptr);
		const auto restore = qScopeGuard([&] { setDefaultEphemerisSource(previous); });
		auto* mov = session.mov();
		QSignalSpy failed(mov, &MovViewModel::runFailed);
		const QVariantMap options{{"startUtc", start}, {"endUtc", end},   {"sampling", 2},
		                          {"rangeNodes", 16},  {"rateNodes", 24}, {"workers", 1},
		                          {"target", 0},       {"boundary", 0},   {"observationSelection", 1}};
		mov->run(options);
		QCOMPARE(page->property("currentView").toInt(), 1);
		QTRY_VERIFY_WITH_TIMEOUT(!mov->busy(), 30000);
		QCOMPARE(failed.count(), 1);
		QCOMPARE(page->property("currentView").toInt(), 0);
		QVERIFY(!mov->error().isEmpty());
		QVERIFY(mov->runs().isEmpty());
		QCOMPARE(form->property("startUtc").toString(), start);
		QCOMPARE(form->property("endUtc").toString(), end);
		QCOMPARE(form->property("sampling").toInt(), 2);
		QCOMPARE(form->property("rangeNodes").toInt(), 16);
		QCOMPARE(form->property("rateNodes").toInt(), 24);
		mov->run(options);
		QVERIFY(session.solution()->publish(source, m_settings.path()));
		QTRY_VERIFY_WITH_TIMEOUT(!mov->busy(), 30000);
		QCOMPARE(failed.count(), 1);
		QVERIFY(mov->error().isEmpty());
		window->close();
	}

	void movLiveWithoutNominal()
	{
		const auto file = qEnvironmentVariable("ORBITAL_UI_MOV_INPUT");
		if (file.isEmpty())
			QSKIP("Live MOV is opt-in with an observation file");
		const auto oldDirectory = getExecutableDir();
		auto* oldSource = getDefaultEphemerisSource();
		auto* oldAsteroids = getDefaultAsteroidSource();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source eph;
		SpkReader asteroids;
		QVERIFY(eph.isAvailable());
		QVERIFY(asteroids.isAvailable());
		setDefaultEphemerisSource(&eph);
		setDefaultAsteroidSource(&asteroids);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(oldSource);
			setDefaultAsteroidSource(oldAsteroids);
			getExecutableDir() = oldDirectory;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		settings.perturbations()->setProperty("enabled", true);
		settings.perturbations()->setUseSb441(true);
		settings.weighting()->setProperty("useWeights", true);
		settings.weighting()->setProperty("observerSigmas", false);
		settings.weighting()->setProperty("timingEnabled", false);
		settings.weighting()->setProperty("debiasEnabled", true);
		settings.nonGrav()->setProperty("enabled", false);
		auto* mov = session.mov();
		mov->selectFile(QUrl::fromLocalFile(file));
		QVERIFY(session.busy());
		QTRY_VERIFY_WITH_TIMEOUT(!mov->busy(), 30000);
		QVERIFY2(mov->error().isEmpty(), qPrintable(mov->error()));
		QVERIFY(!session.solution()->available());
		QVERIFY(mov->externalSource());
		QVERIFY(mov->canRun());
		QVariantMap options{{"startUtc", "2014-01-01 08:00:00"},
		                    {"endUtc", "2014-01-03 00:00:00"},
		                    {"rangeNodes", 32},
		                    {"rateNodes", 32},
		                    {"workers", 4},
		                    {"sampling", 0},
		                    {"target", 0},
		                    {"boundary", 0},
		                    {"observationSelection", 1}};
		mov->run(options);
		QVERIFY(mov->busy());
		QVERIFY(session.busy());
		QVERIFY(!session.solution()->observations()->property("editingEnabled").toBool());
		QTRY_VERIFY_WITH_TIMEOUT(!mov->busy(), 120000);
		QVERIFY2(mov->error().isEmpty(), qPrintable(mov->error()));
		QCOMPARE(mov->runs().size(), 1);
		QCOMPARE(mov->result().value("sampling").toString(), QString("Grid"));
		QVERIFY(mov->result().value("hits").toInt() > 0);
		QVERIFY(mov->result().value("probability").toString().toDouble() > .99);
		QVERIFY(QFileInfo::exists(mov->result().value("report").toString()));
		const auto first = mov->result();
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty()) {
			QDir().mkpath(captures);
			QQmlApplicationEngine engine;
			engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
			engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
			engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
			QVERIFY(!engine.rootObjects().isEmpty());
			auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
			QVERIFY(window);
			window->setProperty("currentPage", 4);
			auto* page = findItem(window->contentItem(), "risk-monitor-page");
			QVERIFY(page);
			page->setProperty("showingMov", true);
			for (const auto size : {QSize(1440, 900), QSize(1024, 720)}) {
				window->resize(size);
				for (int view : {0, 1}) {
					page->setProperty("currentView", view);
					QTest::qWait(150);
					QVERIFY(
					    window->grabWindow().save(captures + QString("/mov-%1-%2.png").arg(size.width()).arg(view)));
				}
			}
			window->close();
		}
		mov->run(options);
		mov->cancel();
		QTRY_VERIFY_WITH_TIMEOUT(!mov->busy(), 30000);
		QVERIFY(!session.busy());
		QVERIFY(!session.solution()->available());
		mov->setSelectedRun(0);
		QCOMPARE(mov->result(), first);
		qInfo("MOV archive: %s", qPrintable(first.value("archive").toString()));
	}

	void riskLiveTc3()
	{
		if (!qEnvironmentVariableIsSet("ORBITAL_UI_LIVE_MC"))
			QSKIP("Live MC is opt-in");
		const auto oldDirectory = getExecutableDir();
		auto* oldSource = getDefaultEphemerisSource();
		auto* oldAsteroids = getDefaultAsteroidSource();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source source;
		SpkReader asteroids;
		QVERIFY(source.isAvailable());
		QVERIFY(asteroids.isAvailable());
		setDefaultEphemerisSource(&source);
		setDefaultAsteroidSource(&asteroids);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(oldSource);
			setDefaultAsteroidSource(oldAsteroids);
			getExecutableDir() = oldDirectory;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		settings.perturbations()->setProperty("enabled", true);
		settings.perturbations()->setPlanetMask(2303);
		settings.perturbations()->setUseSb441(true);
		settings.weighting()->setProperty("useWeights", true);
		settings.weighting()->setProperty("observerSigmas", true);
		settings.weighting()->setProperty("timingEnabled", false);
		settings.weighting()->setProperty("debiasEnabled", true);
		settings.weighting()->setProperty("weightProfile", 1);
		settings.sigmaClip()->setProperty("enabled", true);
		settings.sigmaClip()->setProperty("threshold", 3.0);
		settings.nonGrav()->setProperty("enabled", false);
		QVERIFY(session.openInput(
		    QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/../data/examples/NEO/impacts/2008_TC3.txt")));
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
		QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
		QVERIFY(session.solution()->snapshot()->fit.converged);
		const auto baseline = *session.solution()->snapshot();
		auto* monitor = session.risk();
		monitor->setCenterUtc("2008-10-07 03:00:00");
		monitor->setHalfWindowHours(2);
		QJsonObject covarianceRecord;
		for (int method = 0; method < 2; ++method) {
			monitor->setMethod(method);
			monitor->setSamples(method ? "8" : "32");
			QVERIFY2(monitor->canRun(), qPrintable(monitor->sourceError() + monitor->formError()));
			monitor->run();
			QTRY_VERIFY_WITH_TIMEOUT(!monitor->busy(), 180000);
			QVERIFY2(monitor->error().isEmpty(), qPrintable(monitor->error()));
			qInfo("MC %d: %s; %s; archive=%s", method, qPrintable(monitor->summary()), qPrintable(monitor->warning()),
			      qPrintable(monitor->archivePath()));
			QVERIFY(monitor->events()->count() > 0);
			QFile json(monitor->archivePath() + "/result.json");
			QVERIFY(json.open(QIODevice::ReadOnly));
			const auto record = QJsonDocument::fromJson(json.readAll()).object();
			QCOMPARE(record["completed"].toInt(), method ? 8 : 32);
			QCOMPARE(record["entry_clones"].toInt(), method ? 8 : 32);
			QCOMPARE(record["failed_propagation"].toInt(), 0);
			if (method == 0)
				covarianceRecord = record;
			QCOMPARE(session.solution()->snapshot()->fit.keptObsIndices, baseline.fit.keptObsIndices);
			QCOMPARE(session.solution()->snapshot()->fit.fitEpochOrbit->epochMjd(),
			         baseline.fit.fitEpochOrbit->epochMjd());
			QVERIFY(session.solution()->snapshot()->fit.uncertainty.fitParameterCovariance ==
			        baseline.fit.uncertainty.fitParameterCovariance);
		}
		monitor->setMethod(0);
		monitor->setIntervalMode(1);
		monitor->setRangeEndUtc("2008-10-07 05:00:00");
		QVERIFY2(monitor->canRun(), qPrintable(monitor->sourceError() + monitor->formError()));
		monitor->run();
		QTRY_VERIFY_WITH_TIMEOUT(!monitor->busy(), 180000);
		QVERIFY2(monitor->error().isEmpty(), qPrintable(monitor->error()));
		QFile rangeJson(monitor->archivePath() + "/result.json");
		QVERIFY(rangeJson.open(QIODevice::ReadOnly));
		const auto rangeRecord = QJsonDocument::fromJson(rangeJson.readAll()).object();
		QCOMPARE(rangeRecord["completed"].toInt(), 32);
		QCOMPARE(rangeRecord["entry_clones"].toInt(), 32);
		const auto windowClones = covarianceRecord["clones"].toArray();
		const auto rangeClones = rangeRecord["clones"].toArray();
		for (int i = 0; i < 32; ++i) {
			const auto window = windowClones[i].toObject();
			const auto range = rangeClones[i].toObject();
			QCOMPARE(window["state"], range["state"]);
			const auto entry = [](const QJsonObject& clone) { return clone["entry"].toObject()["epoch"].toDouble(); };
			QVERIFY(std::abs(entry(window) - entry(range)) < 1e-7);
		}
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty()) {
			QQmlApplicationEngine engine;
			engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
			engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
			engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
			QVERIFY(!engine.rootObjects().isEmpty());
			auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
			QVERIFY(window);
			window->setProperty("currentPage", 4);
			QTest::qWait(100);
			auto* page = findItem(window->contentItem(), "risk-monitor-page");
			QVERIFY(page);
			page->setProperty("currentView", 1);
			for (const auto* theme : {"graphite", "light"}) {
				session.setTheme(QString::fromLatin1(theme));
				for (const auto size : {QSize(1440, 900), QSize(1024, 720)}) {
					window->resize(size);
					QTest::qWait(120);
					QVERIFY(window->grabWindow().save(captures + "/tc3-live-" + theme + "-" +
					                                  QString::number(size.width()) + ".png"));
				}
			}
		}
	}

	void riskDefaultWindowCentre()
	{
		namespace time = orbital_core::time;
		SolutionViewModel solution;
		RiskMonitorViewModel monitor(solution);
		QVERIFY(monitor.centerUtc().isEmpty());
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(result, "default-window-centre"));
		QVERIFY(!monitor.centerUtc().isEmpty());
		QVERIFY2(monitor.inputReady(), qPrintable(monitor.formError()));
		const auto before = monitor.draft().value();
		QCOMPARE(time::toMjdTdb(before.startMjdUtc, time::Scale::UTC), 60000.0);
		QCOMPARE(before.endMjdUtc - before.startMjdUtc, 2.0);
		const auto initialCenter = monitor.centerUtc();
		monitor.setIntervalMode(1);
		QCOMPARE(monitor.draft()->startMjdUtc, before.startMjdUtc);
		monitor.setIntervalMode(0);
		monitor.setMethod(1);
		QCOMPARE(monitor.centerUtc(), initialCenter);
		QCOMPARE(monitor.draft()->startMjdUtc, before.startMjdUtc);
		for (const auto* center : {"2025-01-01 00:00:00", "invalid", ""}) {
			monitor.setCenterUtc(QString::fromLatin1(center));
			QVERIFY(solution.publish(result, "retained-window-centre"));
			QCOMPARE(monitor.centerUtc(), QString::fromLatin1(center));
		}
		SolutionViewModel editedSolution;
		RiskMonitorViewModel editedMonitor(editedSolution);
		editedMonitor.setCenterUtc("2025-02-01 12:00:00");
		QVERIFY(editedSolution.publish(result, "pre-edited-window-centre"));
		QCOMPARE(editedMonitor.centerUtc(), QStringLiteral("2025-02-01 12:00:00"));
	}

	void riskIntervalModes()
	{
		namespace time = orbital_core::time;
		SolutionViewModel solution;
		RiskMonitorViewModel monitor(solution);
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(result, "intervals"));
		QCOMPARE(monitor.intervalMode(), 0);
		QVERIFY(monitor.inputReady());
		monitor.setIntervalMode(1);
		QVERIFY2(monitor.inputReady(), qPrintable(monitor.formError()));
		const auto range = monitor.draft().value();
		QCOMPARE(range.intervalMode, RiskMonitorDraft::IntervalMode::DateRange);
		QCOMPARE(time::toMjdTdb(range.startMjdUtc, time::Scale::UTC), 60000.0);
		QCOMPARE(time::toMjdTdb(range.endMjdUtc, time::Scale::UTC), 60000.0 + 30.0 * 365.25);
		QCOMPARE(monitor.rangeYears(), 30.0);
		monitor.setRangeYears(0.5);
		QCOMPARE(monitor.rangeYears(), 0.5);
		const auto end = monitor.rangeEndUtc();
		monitor.setCenterUtc("2023-02-27 00:00:00");
		monitor.setHalfWindowHours(0.5);
		monitor.setIntervalMode(0);
		QVERIFY(monitor.inputReady());
		QCOMPARE(monitor.draft()->intervalMode, RiskMonitorDraft::IntervalMode::EventWindow);
		QVERIFY(std::abs(monitor.draft()->endMjdUtc - monitor.draft()->startMjdUtc - 1.0 / 24.0) < 1e-10);
		monitor.setRangeEndUtc("invalid");
		QVERIFY(monitor.inputReady());
		monitor.setIntervalMode(1);
		QVERIFY(!monitor.inputReady());
		monitor.setRangeEndUtc(end);
		monitor.setCenterUtc("invalid");
		monitor.setHalfWindowHours(0);
		QVERIFY(monitor.inputReady());
		monitor.setIntervalMode(-1);
		monitor.setIntervalMode(2);
		QCOMPARE(monitor.intervalMode(), 1);
		const auto start = monitor.rangeStartUtc();
		monitor.setRangeEndUtc(start);
		QVERIFY(!monitor.inputReady());
		monitor.setRangeEndUtc("2023-02-24 00:00:00");
		QVERIFY(!monitor.inputReady());
		monitor.setRangeEndUtc(end);
		monitor.setRangeStartUtc("2023-02-24 00:00:00");
		QVERIFY(!monitor.inputReady());
		monitor.setRangeStartUtc("2023-02-30 00:00:00");
		QVERIFY(!monitor.inputReady());
		monitor.setRangeStartUtc(start);
		for (double years : {0.0, -1.0, std::numeric_limits<double>::infinity(),
		                     std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::max(), 1e-20}) {
			monitor.setRangeYears(years);
			QVERIFY(!monitor.inputReady());
		}
		monitor.setRangeYears(30);
		QVERIFY(monitor.inputReady());
		monitor.setRangeStartUtc("2023-03-01 01:02:03.123456789");
		monitor.setRangeYears(2);
		QCOMPARE(monitor.rangeYears(), 2.0);
		QCOMPARE(monitor.rangeStartUtc(), QStringLiteral("2023-03-01 01:02:03.123456789"));
		monitor.setRangeEndUtc("2024-03-01 01:02:03.123456789");
		QVERIFY(std::abs(monitor.rangeYears() - 366.0 / 365.25) < 1e-12);
		const auto editedStart = monitor.rangeStartUtc();
		const auto editedEnd = monitor.rangeEndUtc();
		QVERIFY(solution.publish(result, "new-source"));
		QCOMPARE(monitor.rangeStartUtc(), editedStart);
		QCOMPARE(monitor.rangeEndUtc(), editedEnd);
		QVERIFY(monitor.draft()->sourceRevision > range.sourceRevision);
		for (const double epoch : {54745.75, 57753.123456789, 59166.996581801, 60000.999999999}) {
			SolutionViewModel epochSolution;
			RiskMonitorViewModel epochMonitor(epochSolution);
			result.fit.fitEpochOrbit->setEpochMjd(epoch);
			QVERIFY(epochSolution.publish(result, "fractional-epoch"));
			QVERIFY2(epochMonitor.inputReady(), qPrintable(epochMonitor.formError()));
			QCOMPARE(time::toMjdTdb(epochMonitor.draft()->startMjdUtc, time::Scale::UTC), epoch);
			epochMonitor.setIntervalMode(1);
			QVERIFY2(epochMonitor.inputReady(), qPrintable(epochMonitor.formError()));
			const auto interval = epochMonitor.draft().value();
			QCOMPARE(time::toMjdTdb(interval.startMjdUtc, time::Scale::UTC), epoch);
			QCOMPARE(time::toMjdTdb(interval.endMjdUtc, time::Scale::UTC), epoch + 30 * 365.25);
		}
	}

	void ephemeridesRequireConvergedFit()
	{
		EphemerisViewModel ephemerides;
		auto result = fixture();
		ephemerides.setSolution(result, "Fixture", std::nullopt);
		QVERIFY(ephemerides.canGenerate());

		result.fit.converged = false;
		ephemerides.setSolution(result, "Fixture", std::nullopt);
		QVERIFY(ephemerides.available());
		QVERIFY(!ephemerides.canGenerate());
		QCOMPARE(ephemerides.formError(), QStringLiteral("A converged solution is required."));
		ephemerides.generate();
		QCOMPARE(ephemerides.error(), ephemerides.formError());

		result.fit.converged = true;
		ephemerides.setSolution(result, "Fixture", std::nullopt);
		QVERIFY(ephemerides.canGenerate());
	}

	void historicalFitUsesOneUtcDefaultAcrossCalculationForms()
	{
		namespace time = orbital_core::time;
		auto result = fixture();
		result.effectiveConfig.nbody.enabled = true;
		result.fit.fitEpochOrbit->setEpochMjd(22108.546444633510);
		result.fit.orbit->setEpochMjd(52706.0);
		result.observations.objectName = "Historical comet";

		SolutionViewModel solution;
		RiskMonitorViewModel risk(solution);
		EphemerisViewModel ephemerides;
		CloseApproachesViewModel approaches(solution);
		QVERIFY(solution.publish(result, "historical-comet"));
		ephemerides.setSolution(result, "historical-comet", std::nullopt);
		QVERIFY(solution.fitEpochUtc().isEmpty());
		const auto start = solution.calculationStartUtc();
		QVERIFY(!start.isEmpty());
		QCOMPARE(ephemerides.startDate(), start.left(10));
		QVERIFY(ephemerides.startValid());
		QVERIFY(ephemerides.endValid());
		QVERIFY2(risk.inputReady(), qPrintable(risk.formError()));
		QCOMPARE(risk.rangeStartUtc().left(19), start.left(19));

		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import TrajectaUI
            SearchIntervalForm { sourceIdentity: "Historical comet" }
        )",
		                  QUrl());
		std::unique_ptr<QObject> interval(component.create());
		QVERIFY2(interval, qPrintable(component.errorString()));
		QVERIFY(interval->setProperty("seedUtc", start));
		QCOMPARE(interval->property("startUtc").toString(), start);
		const auto end = interval->property("endUtc").toString();
		QVERIFY2(approaches.validationError(start, end, 0.0).isEmpty(),
		         qPrintable(approaches.validationError(start, end, 0.0)));
		QCOMPARE(time::toMjdTdb(risk.draft()->startMjdUtc, time::Scale::UTC), 52706.0);

		result.fit.orbit->setEpochMjd(22109.0);
		result.observations.objectName = "Older comet";
		QVERIFY(solution.publish(result, "older-comet"));
		ephemerides.setSolution(result, "older-comet", std::nullopt);
		const auto today = QDateTime::currentDateTimeUtc().date();
		QCOMPARE(QDate::fromString(ephemerides.startDate(), Qt::ISODate), today);
		QCOMPARE(solution.calculationStartUtc().left(10), today.toString(Qt::ISODate));
		QVERIFY2(risk.inputReady(), qPrintable(risk.formError()));
		QVERIFY(!risk.rangeStartUtc().isEmpty());
	}

	void riskDatesResetOnlyWhenObjectChanges()
	{
		SolutionViewModel solution;
		RiskMonitorViewModel monitor(solution);
		auto first = fixture();
		first.observations.objectName = "First";
		first.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(first, "first-object"));
		monitor.setCenterUtc("2025-02-01 12:00:00");
		monitor.setRangeStartUtc("2025-02-01 00:00:00");
		QVERIFY(solution.publish(first, "first-object-refit"));
		QCOMPARE(monitor.centerUtc(), QStringLiteral("2025-02-01 12:00:00"));
		QCOMPARE(monitor.rangeStartUtc(), QStringLiteral("2025-02-01 00:00:00"));
		auto second = fixture();
		second.observations.objectName = "Second";
		second.effectiveConfig.nbody.enabled = true;
		QVERIFY(solution.publish(second, "second-object"));
		QVERIFY(monitor.centerUtc() != QStringLiteral("2025-02-01 12:00:00"));
		QVERIFY(monitor.rangeStartUtc() != QStringLiteral("2025-02-01 00:00:00"));
		QVERIFY(!monitor.centerUtc().isEmpty());
		QVERIFY(!monitor.rangeStartUtc().isEmpty());
	}

	void approachIntervalResetsOnlyWhenObjectChanges()
	{
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import TrajectaUI
            SearchIntervalForm { seedUtc: "2023-02-24 00:00:00"; sourceIdentity: "First" }
        )",
		                  QUrl());
		std::unique_ptr<QObject> form(component.create());
		QVERIFY2(form, qPrintable(component.errorString()));
		QCOMPARE(form->property("startUtc").toString(), QStringLiteral("2023-02-24 00:00:00"));
		QVERIFY(form->setProperty("startUtc", "2025-01-01 00:00:00"));
		QVERIFY(form->setProperty("seedUtc", "2024-02-24 00:00:00"));
		QCOMPARE(form->property("startUtc").toString(), QStringLiteral("2025-01-01 00:00:00"));
		QVERIFY(form->setProperty("sourceIdentity", "Second"));
		QTRY_COMPARE(form->property("startUtc").toString(), QStringLiteral("2024-02-24 00:00:00"));
	}

	void orbitSceneSnapshotAndTime()
	{
		OrbitSceneModel model;
		QVERIFY(!model.available());
		const auto result = fixture();
		const auto orbit = result.fit.orbit;
		model.setOrbit(orbit);
		QVERIFY(model.available());
		QVERIFY(model.canAnimate());
		QCOMPARE(model.epochMjdTdb(), 60000.0);
		QCOMPARE(model.epochUtc().left(10), QStringLiteral("2023-02-24"));
		const auto utcTime = QTime::fromString(model.epochUtc().mid(11), QStringLiteral("HH:mm:ss.zzz"));
		QVERIFY(utcTime.isValid());
		QVERIFY(std::abs(utcTime.msecsTo(QTime(23, 58, 50, 816))) <= 3);
		QCOMPARE(model.spanDays(), orbit->getOrbitalPeriod());
		const auto original = model.objectPosition();
		model.setPosition(0.25);
		QVERIFY((original - model.objectPosition()).length() > 0.1);
		QCOMPARE(model.epochMjdTdb(), 60000 + orbit->getOrbitalPeriod() * 0.25);
		QVERIFY(model.epochUtc().left(10) != QStringLiteral("2023-02-24"));
		model.setPosition(1);
		QVERIFY((original - model.objectPosition()).length() < 1e-5);
		model.setPosition(2);
		QCOMPARE(model.position(), 1.0);
		model.setPosition(std::numeric_limits<double>::quiet_NaN());
		QCOMPARE(model.position(), 1.0);
		model.setPosition(-1);
		QCOMPARE(model.objectPosition(), original);
		QCOMPARE(orbit->epochMjd(), 60000.0);
		for (double e : {1.0, 1.2}) {
			orbital_core::cometary::Vector6 elements;
			elements << 0.8, e, 0.4, 1.2, 0.8, 0.2;
			const auto state = orbital_core::cometary::toCartesian(elements).value();
			model.setOrbit(
			    make_orbit_from_states_vectors(StateVectors(Vec3D(state(0), state(1), state(2)),
			                                                Vec3D(state(3), state(4), state(5)) / Constants::GAUSS_K),
			                                   2460000.5));
			QCOMPARE(model.canAnimate(), e > 1.0);
			model.setPosition(0.5);
			QVERIFY(std::isfinite(model.objectPosition().length()));
			if (e == 1.0)
				QCOMPARE(model.position(), 0.0);
		}
		auto historical = std::make_shared<Orbit>(*orbit);
		historical->setEpochMjd(22108.546444633510);
		model.setOrbit(historical);
		QVERIFY(model.epochUtc().isEmpty());
		model.setOrbit(nullptr);
		QVERIFY(model.epochUtc().isEmpty());
		QVERIFY(!model.available());
		QVERIFY(!model.canAnimate());
		QVERIFY(model.planets().isEmpty());
	}

	void nativeOrbitScene()
	{
		if (!qEnvironmentVariableIsSet("ORBITAL_UI_3D"))
			QSKIP("Native 3D rendering is opt-in");
		const auto oldExecutableDir = getExecutableDir();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source de440;
		auto* oldSource = getDefaultEphemerisSource();
		setDefaultEphemerisSource(&de440);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(oldSource);
			getExecutableDir() = oldExecutableDir;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto result = fixture();
		QVERIFY(session.solution()->publish(result, "test"));
		auto* sceneModel = session.solution()->scene();
		QCOMPARE(sceneModel->planets().size(), 8);
		QVERIFY(sceneModel->error().isEmpty());
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QVERIFY(QTest::qWaitForWindowExposed(window));
		const auto click = [&](const QString& name) {
			auto* item = findItem(window->contentItem(), name);
			QVERIFY2(item, qPrintable(name));
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
			                  item->mapToScene({item->width() / 2, item->height() / 2}).toPoint());
		};
		click("navigation-page-2");
		QTRY_VERIFY(findItem(window->contentItem(), "orbit-scene"));
		auto* scene = findItem(window->contentItem(), "orbit-scene");
		QVERIFY(!scene->findChild<QObject*>("orbit-body-Pluto"));
		QVERIFY(!scene->findChild<QObject*>("orbit-path-Pluto"));
		auto* earth = scene->findChild<QObject*>("orbit-body-Earth");
		QVERIFY(earth);
		QQmlListReference earthMaterials(earth, "materials");
		QCOMPARE(earthMaterials.count(), 1);
		auto* earthMaterial = earthMaterials.at(0);
		QVERIFY(!earthMaterial->property("baseColorMap").value<QObject*>());
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			QTRY_COMPARE(earthMaterial->property("baseColor").value<QColor>(),
			             QColor(theme == "light" ? "#4877b9" : "#6f9fdc"));
			for (const QSize size : {QSize(1024, 720), QSize(1512, 600), QSize(1512, 982), QSize(2554, 1312)}) {
				window->resize(size);
				QTest::qWait(250);
				const auto bounds = scene->mapRectToScene(scene->boundingRect());
				QVERIFY(bounds.width() > 800);
				QVERIFY(bounds.height() > 450);
				QVERIFY(bounds.right() <= window->width());
				auto* legend = findItem(scene, "orbit-legend");
				QVERIFY(legend);
				const auto legendBounds = legend->mapRectToScene(legend->boundingRect());
				QVERIFY(bounds.contains(legendBounds));
				for (auto* entry : legend->childItems()) {
					if (entry->objectName().startsWith("orbit-follow-")) {
						QVERIFY(legendBounds.contains(entry->mapRectToScene(entry->boundingRect())));
						QVERIFY(entry->width() > 35);
						QCOMPARE(entry->height(), 28.0);
					}
				}
				auto* playback = findItem(window->contentItem(), "orbit-time");
				QVERIFY(playback && playback->width() > 200);
				auto* utc = findItem(window->contentItem(), "orbit-utc");
				QVERIFY(utc && utc->isVisible());
				QCOMPARE(utc->property("value").toString(), sceneModel->epochUtc());
				QVERIFY(utc->mapRectToScene(utc->boundingRect()).right() <=
				        playback->mapRectToScene(playback->boundingRect()).left());
				const auto image = window->grabWindow();
				QVERIFY(!image.isNull());
				const auto dpr = image.devicePixelRatio();
				const auto area = QRect(qRound((bounds.x() + 30) * dpr), qRound((bounds.y() + 160) * dpr),
				                        qRound((bounds.width() - 60) * dpr), qRound((bounds.height() - 240) * dpr));
				int colored = 0;
				for (int y = area.top(); y < area.bottom(); y += 2) {
					for (int x = area.left(); x < area.right(); x += 2) {
						const auto c = image.pixelColor(x, y);
						if (std::abs(c.green() - c.red()) > 35 || std::abs(c.blue() - c.red()) > 35)
							++colored;
					}
				}
				QVERIFY2(colored > 120, qPrintable(QString::number(colored)));
				if (!captures.isEmpty())
					QVERIFY(image.save(captures + "/orbit-" + theme + "-" + QString::number(size.width()) + ".png"));
			}
		}
		window->resize(1512, 982);
		session.setTheme("graphite");
		QTest::qWait(120);
		const double distance = scene->property("cameraDistance").toDouble();
		click("orbit-zoom-in");
		QVERIFY(scene->property("cameraDistance").toDouble() < distance);
		click("orbit-reset-view");
		QCOMPARE(scene->property("cameraDistance").toDouble(), distance);
		const auto before = window->grabWindow();
		click("orbit-top-view");
		QTest::qWait(100);
		QCOMPARE(scene->property("cameraRotation").value<QVector3D>().x(), -90.0f);
		QVERIFY(window->grabWindow() != before);
		const auto center = scene->mapToScene({scene->width() * 0.5, scene->height() * 0.5}).toPoint();
		QTest::mousePress(window, Qt::LeftButton, Qt::NoModifier, center);
		QTest::mouseMove(window, center + QPoint(100, 50), 40);
		QTest::qWait(80);
		QTest::mouseMove(window, center + QPoint(180, 100), 40);
		QTest::qWait(80);
		QTest::mouseRelease(window, Qt::LeftButton, Qt::NoModifier, center + QPoint(180, 100));
		QVERIFY(scene->property("cameraRotation").value<QVector3D>() != QVector3D(-90, 0, 0));
		click("orbit-reset-view");
		const auto initialPosition = sceneModel->objectPosition();
		click("orbit-play");
		QTRY_VERIFY(sceneModel->position() > 0.001);
		click("orbit-play");
		QVERIFY((initialPosition - sceneModel->objectPosition()).length() > 0.001);
		click("orbit-reset-time");
		QCOMPARE(sceneModel->position(), 0.0);
		click("orbit-time");
		QVERIFY(sceneModel->position() > 0.45 && sceneModel->position() < 0.55);
		click("orbit-reset-time");
		const auto planets = sceneModel->planets();
		for (qsizetype index = 0; index < planets.size(); ++index) {
			const auto name = planets[index].toMap().value("name").toString();
			auto* path = scene->findChild<QObject*>("orbit-path-" + name);
			QVERIFY(path);
			QVERIFY(path->property("color").value<QColor>().lightnessF() > 0.15);
			click("orbit-follow-" + name);
			QCOMPARE(scene->property("selectedTarget").toString(), name);
			QTRY_VERIFY((scene->property("cameraCenter").value<QVector3D>() -
			             sceneModel->planetPositions()[index].value<QVector3D>())
			                .length() < 1e-5f);
			QTRY_VERIFY(std::abs(scene->property("cameraDistance").toDouble() - 1.0) < 1e-5);
			QVERIFY(findItem(scene, "orbit-follow-" + name)->property("checked").toBool());
			sceneModel->setPosition(0.03);
			QTRY_VERIFY((scene->property("cameraCenter").value<QVector3D>() -
			             sceneModel->planetPositions()[index].value<QVector3D>())
			                .length() < 1e-5f);
			sceneModel->setPosition(0);
		}
		click("orbit-follow-object");
		QTRY_VERIFY((scene->property("cameraCenter").value<QVector3D>() - sceneModel->objectPosition()).length() <
		            1e-5f);
		click("orbit-play");
		QTRY_VERIFY(sceneModel->position() > 0.002);
		QVERIFY((scene->property("cameraCenter").value<QVector3D>() - sceneModel->objectPosition()).length() < 1e-5f);
		click("orbit-play");
		click("orbit-top-view");
		QCOMPARE(scene->property("selectedTarget").toString(), QString("object"));
		QVERIFY((scene->property("cameraCenter").value<QVector3D>() - sceneModel->objectPosition()).length() < 1e-5f);
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/orbit-follow-object.png"));
		auto* cameraOrigin = scene->findChild<QObject*>("orbit-camera-origin");
		QVERIFY(cameraOrigin);
		const QVector3D panOffset(0.1f, -0.2f, 0.05f);
		QVERIFY(cameraOrigin->setProperty("position", panOffset));
		sceneModel->setPosition(0.1);
		QVERIFY(
		    (scene->property("cameraCenter").value<QVector3D>() - sceneModel->objectPosition() - panOffset).length() <
		    1e-5f);
		click("orbit-follow-object");
		QTRY_VERIFY((scene->property("cameraCenter").value<QVector3D>() - sceneModel->objectPosition()).length() <
		            1e-5f);
		click("orbit-follow-perihelion");
		QTest::qWait(550);
		const auto perihelionCenter = scene->property("cameraCenter").value<QVector3D>();
		QVERIFY(perihelionCenter.length() > 0);
		sceneModel->setPosition(0.2);
		QCOMPARE(scene->property("cameraCenter").value<QVector3D>(), perihelionCenter);
		click("orbit-reset-view");
		QCOMPARE(scene->property("selectedTarget").toString(), QString("sun"));
		QCOMPARE(scene->property("cameraCenter").value<QVector3D>(), QVector3D());
		QCOMPARE(scene->property("cameraDistance").toDouble(), distance);
		click("orbit-reset-time");
		for (const QString theme : {"graphite", "observatory", "light"}) {
			session.setTheme(theme);
			QTest::qWait(100);
			if (!captures.isEmpty())
				QVERIFY(window->grabWindow().save(captures + "/orbit-theme-" + theme + ".png"));
		}
		session.setTheme("graphite");
		QTest::qWait(100);
		auto* viewport = findItem(scene, "orbit-viewport");
		auto* gridViewport = findItem(scene, "orbit-grid-viewport");
		QVERIFY(viewport && gridViewport);
		const double initialGridStep = scene->property("gridStep").toDouble();
		QVERIFY(QMetaObject::invokeMethod(scene, "zoom", Q_ARG(QVariant, 25.0)));
		QVERIFY(scene->property("gridStep").toDouble() > initialGridStep);
		QVERIFY(gridViewport->property("halfDivisions").toInt() * gridViewport->property("step").toDouble() >= 50.0);
		QTest::qWait(100);
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/orbit-grid-zoomed-out.png"));
		QVERIFY(QMetaObject::invokeMethod(scene, "resetView"));
		QCOMPARE(scene->property("gridStep").toDouble(), initialGridStep);
		for (const float elevation : {-15.0f, 15.0f}) {
			QVERIFY(cameraOrigin->setProperty("eulerRotation", QVariant::fromValue(QVector3D(elevation, 0, 0))));
			QTest::qWait(100);
			const auto compositeWithGrid = window->grabWindow();
			if (!captures.isEmpty())
				QVERIFY(compositeWithGrid.save(captures + "/orbit-grid-" + QString::number(elevation) + ".png"));
			const auto withGrid = viewport->grabToImage();
			QVERIFY(withGrid);
			QTRY_VERIFY(!withGrid->image().isNull());
			QVERIFY(gridViewport->setProperty("visible", false));
			QTest::qWait(100);
			QVERIFY(window->grabWindow() != compositeWithGrid);
			const auto withoutGrid = viewport->grabToImage();
			QVERIFY(withoutGrid);
			QTRY_VERIFY(!withoutGrid->image().isNull());
			const auto foreground = withGrid->image();
			QCOMPARE(foreground, withoutGrid->image());
			int opaquePixels = 0;
			int transparentPixels = 0;
			for (int y = 0; y < foreground.height(); ++y) {
				for (int x = 0; x < foreground.width(); ++x) {
					const auto alpha = qAlpha(foreground.pixel(x, y));
					if (alpha == 255)
						++opaquePixels;
					if (alpha == 0)
						++transparentPixels;
				}
			}
			QVERIFY(opaquePixels > 100);
			QVERIFY(transparentPixels > foreground.width() * foreground.height() / 2);
			QVERIFY(gridViewport->setProperty("visible", true));
		}
		click("orbit-reset-view");
		click("orbit-layers");
		QTest::qWait(100);
		click("orbit-layer-grid");
		QVERIFY(!scene->property("showGrid").toBool());
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/orbit-layers.png"));
		QTest::keyClick(window, Qt::Key_Escape);
		QTest::qWait(100);
		click("orbit-play");
		QTRY_VERIFY(sceneModel->position() > 0);
		click("navigation-page-0");
		const double stopped = sceneModel->position();
		QTest::qWait(150);
		QCOMPARE(sceneModel->position(), stopped);
		QVERIFY(!findItem(window->contentItem(), "orbit-scene"));
		window->close();
	}

	void bundledFontsAvailable()
	{
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import TrajectaUI
            import "shared/UI.js" as UI
            QtObject {
                property string textFamily: AppTheme.textFont.name
                property string italicFamily: AppTheme.textItalicFont.name
                property string semiBoldFamily: AppTheme.textSemiBoldFont.name
                property string boldFamily: AppTheme.textBoldFont.name
                property string monoFamily: AppTheme.monoFont.name
                property string monoBoldFamily: AppTheme.monoBoldFont.name
                property string uiTextFamily: UI.font
                property string uiMonoFamily: UI.monoFont
            }
        )",
		                  QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/font-check.qml"));
		std::unique_ptr<QObject> root(component.create());
		QVERIFY2(root, qPrintable(component.errorString()));
#ifdef Q_OS_WIN
		QCOMPARE(root->property("textFamily").toString(), QStringLiteral("Inter"));
		QCOMPARE(root->property("monoFamily").toString(), QStringLiteral("JetBrains Mono"));
		QCOMPARE(root->property("monoBoldFamily").toString(), QStringLiteral("JetBrains Mono"));
		const QString textFamily = root->property("uiTextFamily").toString();
		QVERIFY(textFamily == QStringLiteral("Segoe UI Variable Text") || textFamily == QStringLiteral("Segoe UI") ||
		        textFamily == root->property("textFamily").toString());
		const QString monoFamily = root->property("uiMonoFamily").toString();
		QVERIFY(monoFamily == QStringLiteral("Cascadia Mono") || monoFamily == QStringLiteral("Consolas") ||
		        monoFamily == root->property("monoFamily").toString());
#else
		QCOMPARE(root->property("textFamily").toString(), QStringLiteral("Source Sans 3"));
		QCOMPARE(root->property("italicFamily").toString(), QStringLiteral("Source Sans 3"));
		QCOMPARE(root->property("semiBoldFamily").toString(), QStringLiteral("Source Sans 3"));
		QCOMPARE(root->property("boldFamily").toString(), QStringLiteral("Source Sans 3"));
		QCOMPARE(root->property("monoFamily").toString(), QStringLiteral("Cascadia Mono"));
		QCOMPARE(root->property("uiTextFamily").toString(), root->property("textFamily").toString());
		QCOMPARE(root->property("uiMonoFamily").toString(), root->property("monoFamily").toString());
#endif
	}

	void bundledApplicationIcon()
	{
		const QString path = QStringLiteral(":/qt/qml/TrajectaUI/resources/app.ico");
		QImageReader reader(path);
		QCOMPARE(reader.imageCount(), 9);
		const QIcon icon(path);
		QVERIFY(!icon.isNull());
		const QIcon windowIcon(QStringLiteral(":/qt/qml/TrajectaUI/resources/icons/workspace/app-logo.svg"));
		QVERIFY(!windowIcon.pixmap(32, 32).isNull());
		for (int size : {16, 20, 24, 32, 40, 48, 64, 128, 256}) {
			const auto image = icon.pixmap(QSize(size, size), 1.0).toImage();
			QCOMPARE(image.size(), QSize(size, size));
			QCOMPARE(image.pixelColor(0, 0).alpha(), 0);
			int colored = 0;
			for (int y = 0; y < size; ++y) {
				for (int x = 0; x < size; ++x) {
					const auto pixel = image.pixelColor(x, y);
					if (pixel.alpha() > 200 && pixel.green() > pixel.red() && pixel.blue() > pixel.red())
						++colored;
				}
			}
			QVERIFY(colored > size * size / 30);
		}
	}

	void modalBackdropState()
	{
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData(R"(
            import QtQuick
            import QtQuick.Controls.Basic
            import TrajectaUI
            ApplicationWindow {
                id: host
                width: 600; height: 500; visible: true
                readonly property bool overlayVisible: Overlay.overlay.visible
                DialogBackdrop { objectName: "backdrop"; anchors.fill: parent; active: DialogState.activeFor(host) }
                ToolTip { objectName: "hint"; text: "Context help"; timeout: -1 }
                BaseDialog { objectName: "first"; title: "First" }
                BaseDialog { objectName: "second"; title: "Second" }
            }
        )",
		                  QUrl());
		std::unique_ptr<QObject> root(component.create());
		QVERIFY2(root, qPrintable(component.errorString()));
		std::unique_ptr<QObject> other(component.create());
		QVERIFY2(other, qPrintable(component.errorString()));
		auto* backdrop = root->findChild<QObject*>("backdrop");
		auto* otherBackdrop = other->findChild<QObject*>("backdrop");
		auto* hint = root->findChild<QObject*>("hint");
		auto* first = root->findChild<QObject*>("first");
		auto* second = root->findChild<QObject*>("second");
		QVERIFY(backdrop && otherBackdrop && hint && first && second);
		QVERIFY(!backdrop->property("active").toBool());
		QVERIFY(QMetaObject::invokeMethod(hint, "open"));
		QTRY_VERIFY(root->property("overlayVisible").toBool());
		QVERIFY(!backdrop->property("active").toBool());
		QVERIFY(QMetaObject::invokeMethod(hint, "close"));
		QVERIFY(QMetaObject::invokeMethod(first, "open"));
		QTRY_VERIFY(backdrop->property("active").toBool());
		QVERIFY(!otherBackdrop->property("active").toBool());
		QVERIFY(QMetaObject::invokeMethod(second, "open"));
		QVERIFY(QMetaObject::invokeMethod(first, "close"));
		QTRY_VERIFY(!first->property("visible").toBool());
		QVERIFY(backdrop->property("active").toBool());
		QVERIFY(QMetaObject::invokeMethod(second, "close"));
		QTRY_VERIFY(!backdrop->property("active").toBool());
		first->setProperty("modal", false);
		QVERIFY(QMetaObject::invokeMethod(first, "open"));
		QTRY_VERIFY(first->property("visible").toBool());
		QVERIFY(!backdrop->property("active").toBool());
		first->setProperty("modal", true);
		QTRY_VERIFY(backdrop->property("active").toBool());
		delete first;
		QTRY_VERIFY(!backdrop->property("active").toBool());
	}

	void snapshotAndResiduals()
	{
		SolutionViewModel model;
		QVERIFY(model.publish(fixture(), "saved"));
		QCOMPARE(model.total(), 3);
		QCOMPARE(model.included(), 2);
		QCOMPARE(model.arcDays(), 2.0);
		QCOMPARE(model.stationCount(), 1);
		QCOMPARE(model.epochMjd(), QString("60000.000000000"));
		QCOMPARE(model.stateVector().size(), 6);
		auto physicalState = fixture();
		physicalState.fit.orbit->setV({1, 2, 3});
		model.publish(physicalState, "saved");
		QCOMPARE(model.stateVector()[3].toDouble(), Constants::GAUSS_K);
		QScatterSeries series;
		model.observations()->fillSeries(&series, true, false, false);
		QCOMPARE(series.count(), 2);
		orbital_core::FitPipelineResult failed;
		failed.fit.cancelled = true;
		QVERIFY(!model.publish(failed, "failed"));
		QCOMPARE(model.objectTitle(), QString("Fixture"));
		QCOMPARE(model.directory(), QString("saved"));
		QCOMPARE(model.observations()->rowCount(), 3);
	}

	void cancellationBeforeWorkerStarts()
	{
		OrbitWorker worker;
		QSignalSpy errors(&worker, &OrbitWorker::errorOccurred);
		QSignalSpy completed(&worker, &OrbitWorker::optimizationFinished);
		std::stop_source cancellation;
		cancellation.request_stop();
		worker.process({{}, QString("missing.xml")}, cancellation.get_token());
		QCOMPARE(errors.count(), 1);
		QCOMPARE(errors.front().front().toString(), QString("Cancelled"));
		QCOMPARE(completed.count(), 0);
	}

	void manualObservationEditsAreSeparateFromAcceptedResults()
	{
		auto result = fixture();
		ObservationTableModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		model.setResult(result);
		model.sortColumn(1, true);
		model.setStationFilter("A01");
		model.selectRow(0);
		QCOMPARE(model.selectedSourceRow(), 2);
		model.setIncluded(model.selectedSourceRow(), false);
		QCOMPARE(model.pendingCount(), 1);
		QCOMPARE(model.manualCount(), 1);
		QCOMPARE(model.overrides().front().id, result.observations.all[2].id);
		QVERIFY(model.selection().value("included").toBool());
		QVERIFY(!model.selection().value("desiredIncluded").toBool());
		QCOMPARE(model.stationResiduals()->included(), 2);
		QScatterSeries series;
		model.fillSeries(&series, true, false, false);
		QCOMPARE(series.count(), 2);
		model.setEditingEnabled(false);
		model.setIncluded(0, false);
		model.discardEdits();
		model.resetManualSelection();
		QCOMPARE(model.pendingCount(), 1);
		model.setEditingEnabled(true);
		model.setIncluded(2, true);
		QCOMPARE(model.pendingCount(), 0);
		QCOMPARE(model.manualCount(), 0);
		model.setIncluded(1, true);
		QCOMPARE(model.overrides().front().disposition, orbital_core::ObservationDisposition::Include);
		model.discardEdits();
		QVERIFY(model.overrides().empty());
		result.observationOverrides = {{result.observations.all[1].id, orbital_core::ObservationDisposition::Exclude}};
		model.setResult(result);
		QCOMPARE(model.pendingCount(), 0);
		QCOMPARE(model.manualCount(), 1);
		model.resetManualSelection();
		QCOMPARE(model.pendingCount(), 1);
		QCOMPARE(model.manualCount(), 0);
		model.discardEdits();
		QCOMPARE(model.pendingCount(), 0);
		QCOMPARE(model.manualCount(), 1);
		result.fit.keptObsIndices.clear();
		result.fit.residuals = result.fit.allResiduals;
		result.fit.allResiduals.clear();
		result.observationOverrides.clear();
		model.setResult(result);
		QCOMPARE(model.stationResiduals()->included(), 3);
		model.selectSourceRow(1);
		QCOMPARE(model.selection().value("raResidual").toString(), QString("0.100"));
	}

	void liveBorisovSettingsRoundTrips_data()
	{
		QTest::addColumn<bool>("manualSelection");
		QTest::newRow("automatic") << false;
		QTest::newRow("manual-include-exclude") << true;
	}

	void liveBorisovSettingsRoundTrips()
	{
		if (!qEnvironmentVariableIsSet("ORBITAL_UI_LIVE_FIT"))
			QSKIP("Live fit is opt-in");
		QFETCH(bool, manualSelection);
		const auto oldDirectory = getExecutableDir();
		auto* oldSource = getDefaultEphemerisSource();
		auto* oldAsteroids = getDefaultAsteroidSource();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source source;
		SpkReader asteroids;
		QVERIFY(source.isAvailable());
		QVERIFY(asteroids.isAvailable());
		setDefaultEphemerisSource(&source);
		setDefaultAsteroidSource(&asteroids);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(oldSource);
			setDefaultAsteroidSource(oldAsteroids);
			getExecutableDir() = oldDirectory;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		settings.perturbations()->setProperty("enabled", true);
		settings.perturbations()->setPlanetMask(2303);
		settings.perturbations()->setUseSb441(true);
		settings.perturbations()->setProperty("integratorInitialStepDays", 10.0);
		settings.weighting()->setProperty("useWeights", true);
		settings.weighting()->setProperty("observerSigmas", true);
		settings.weighting()->setProperty("timingEnabled", false);
		settings.weighting()->setProperty("debiasEnabled", true);
		settings.weighting()->setProperty("weightProfile", 1);
		settings.sigmaClip()->setProperty("enabled", true);
		settings.sigmaClip()->setProperty("threshold", 3.0);
		settings.nonGrav()->setProperty("model", 1);
		settings.nonGrav()->setProperty("solveMask", 15);
		settings.nonGrav()->setProperty("enabled", true);
		QTemporaryDir inputs;
		const auto file = inputs.filePath("borisov.xml");
		QVERIFY(
		    QFile::copy(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/../data/examples/comets/2I_Borisov_ADES.xml", file));
		QVERIFY(session.openInput(QUrl::fromLocalFile(file)));
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
		QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
		QVERIFY(session.solution()->available());
		QVERIFY(QFile::remove(file));
		const auto run = [&](const char* stage) {
			session.run();
			QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
			QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
			QVERIFY(session.solution()->snapshot()->fit.converged);
			const auto& fit = session.solution()->snapshot()->fit;
			qInfo("%s: kept=%zu RMS=%.17g A1=%.17g A2=%.17g A3=%.17g DT=%.17g", stage, fit.keptObsIndices.size(),
			      fit.cost, fit.finalNonGrav.A1(), fit.finalNonGrav.A2(), fit.finalNonGrav.A3(), fit.finalNonGrav.DT());
		};
		auto* rows = session.solution()->observations();
		const auto automaticBaseline = *session.solution()->snapshot();
		if (manualSelection) {
			const auto& fit = session.solution()->snapshot()->fit;
			rows->setIncluded(static_cast<int>(fit.keptObsIndices.at(fit.keptObsIndices.size() / 2)), false);
			const auto& residuals = fit.allResiduals.data();
			for (std::size_t i = 0; i < residuals.size(); ++i) {
				const auto& weight = residuals[i].astrometricWeight;
				if (weight && weight->postFitChiSquared && *weight->postFitChiSquared < 12.0 &&
				    std::ranges::find(fit.keptObsIndices, i) == fit.keptObsIndices.end()) {
					rows->setIncluded(static_cast<int>(i), true);
					break;
				}
			}
			QCOMPARE(rows->manualCount(), 2);
			run("manual baseline");
			QVERIFY(!QTest::currentTestFailed());
		}
		const auto baseline = *session.solution()->snapshot();
		const auto verifyReturn = [&](const orbital_core::FitPipelineResult& expected) {
			const auto& actual = *session.solution()->snapshot();
			QVERIFY(actual.observationOverrides == expected.observationOverrides);
			QCOMPARE(rows->manualCount(), static_cast<int>(expected.observationOverrides.size()));
			QCOMPARE(rows->pendingCount(), 0);
			QCOMPARE(actual.observations.all.size(), expected.observations.all.size());
			for (std::size_t i = 0; i < actual.observations.all.size(); ++i)
				QVERIFY(actual.observations.all[i].id == expected.observations.all[i].id);
			const auto& a = actual.fit;
			const auto& b = expected.fit;
			QVERIFY(a.converged && b.converged);
			QVERIFY(a.keptObsIndices == b.keptObsIndices);
			QVERIFY(a.cost == b.cost);
			for (const auto& [left, right] :
			     {std::pair{a.fitEpochOrbit, b.fitEpochOrbit}, std::pair{a.orbit, b.orbit}}) {
				QVERIFY(left->getEpoch() == right->getEpoch());
				QVERIFY(left->getR().x == right->getR().x && left->getR().y == right->getR().y &&
				        left->getR().z == right->getR().z);
				QVERIFY(left->getV().x == right->getV().x && left->getV().y == right->getV().y &&
				        left->getV().z == right->getV().z);
			}
			QVERIFY(a.finalNonGrav.A1() == b.finalNonGrav.A1());
			QVERIFY(a.finalNonGrav.A2() == b.finalNonGrav.A2());
			QVERIFY(a.finalNonGrav.A3() == b.finalNonGrav.A3());
			QVERIFY(a.finalNonGrav.DT() == b.finalNonGrav.DT());
			QVERIFY(a.finalNonGrav.enabled() == b.finalNonGrav.enabled());
			const auto& ac = a.uncertainty;
			const auto& bc = b.uncertainty;
			QCOMPARE(ac.fitParameterCovariance.rows(), bc.fitParameterCovariance.rows());
			QVERIFY((ac.fitParameterCovariance.array() == bc.fitParameterCovariance.array()).all());
			QVERIFY((ac.fitParameterCovarianceFormal.array() == bc.fitParameterCovarianceFormal.array()).all());
			QVERIFY(ac.covarianceObjectiveScale == bc.covarianceObjectiveScale);
			QVERIFY(ac.cometaryElementCovariance && bc.cometaryElementCovariance);
			QVERIFY(
			    (ac.cometaryElementCovariance->covariance.array() == bc.cometaryElementCovariance->covariance.array())
			        .all());
			QVERIFY((ac.cometaryElementCovariance->formalCovariance.array() ==
			         bc.cometaryElementCovariance->formalCovariance.array())
			            .all());
			QCOMPARE(a.allResiduals.count(), b.allResiduals.count());
			for (std::size_t i = 0; i < a.allResiduals.count(); ++i) {
				const auto& ar = a.allResiduals.data()[i];
				const auto& br = b.allResiduals.data()[i];
				QVERIFY(ar.dRA_arcsec == br.dRA_arcsec && ar.dDEC_arcsec == br.dDEC_arcsec);
				QVERIFY(ar.astrometricWeight && br.astrometricWeight);
				QVERIFY(ar.astrometricWeight->sigmaRaArcsec == br.astrometricWeight->sigmaRaArcsec);
				QVERIFY(ar.astrometricWeight->sigmaDecArcsec == br.astrometricWeight->sigmaDecArcsec);
				QVERIFY(ar.astrometricWeight->postFitChiSquared == br.astrometricWeight->postFitChiSquared);
			}
			qInfo("Return matches mask, states, non-grav, residuals and full covariance exactly");
		};
		settings.sigmaClip()->setProperty("enabled", false);
		run("clipping disabled");
		QVERIFY(!QTest::currentTestFailed());
		const auto unclipped = *session.solution()->snapshot();
		QCOMPARE(unclipped.fit.keptObsIndices.size(),
		         unclipped.observations.all.size() - static_cast<std::size_t>(manualSelection));
		QCOMPARE(unclipped.fit.sigmaClipPassesRun, 0);
		for (std::size_t i = 0; i < unclipped.observations.all.size(); ++i) {
			const auto id = unclipped.observations.all[i].id;
			const bool excluded = std::ranges::any_of(unclipped.observationOverrides, [id](const auto& entry) {
				return entry.id == id && entry.disposition == orbital_core::ObservationDisposition::Exclude;
			});
			QCOMPARE(std::ranges::find(unclipped.fit.keptObsIndices, i) != unclipped.fit.keptObsIndices.end(),
			         !excluded);
		}
		settings.sigmaClip()->setProperty("enabled", true);
		run("clipping enabled restored");
		verifyReturn(baseline);
		QVERIFY(!QTest::currentTestFailed());
		settings.sigmaClip()->setProperty("enabled", false);
		run("clipping disabled restored");
		verifyReturn(unclipped);
		QVERIFY(!QTest::currentTestFailed());
		settings.sigmaClip()->setProperty("enabled", true);
		settings.sigmaClip()->setProperty("threshold", 2.0);
		run("sigma 2");
		QVERIFY(!QTest::currentTestFailed());
		settings.sigmaClip()->setProperty("threshold", 3.0);
		run("sigma 3 restored");
		verifyReturn(baseline);
		QVERIFY(!QTest::currentTestFailed());
		settings.weighting()->setProperty("weightProfile", 2);
		run("uniform weights");
		QVERIFY(!QTest::currentTestFailed());
		QVERIFY(session.solution()->snapshot()->selection == orbital_core::FitSelectionMode::Configured);
		QVERIFY(session.solution()->snapshot()->fit.sigmaClipPassesRun > 0);
		QVERIFY(session.solution()->snapshot()->fit.converged);
		settings.weighting()->setProperty("weightProfile", 1);
		run("vfcc17 restored");
		verifyReturn(baseline);
		QVERIFY(!QTest::currentTestFailed());
		settings.nonGrav()->setProperty("enabled", false);
		run("non-grav off");
		QVERIFY(!QTest::currentTestFailed());
		const auto gravity = *session.solution()->snapshot();
		QCOMPARE(gravity.fit.uncertainty.fitParameterCovariance.rows(), 6);
		settings.sigmaClip()->setProperty("enabled", false);
		run("gravity clipping disabled");
		QVERIFY(!QTest::currentTestFailed());
		QCOMPARE(session.solution()->snapshot()->fit.keptObsIndices.size(),
		         gravity.observations.all.size() - static_cast<std::size_t>(manualSelection));
		settings.sigmaClip()->setProperty("enabled", true);
		run("gravity clipping enabled restored");
		verifyReturn(gravity);
		QVERIFY(!QTest::currentTestFailed());
		settings.nonGrav()->setProperty("enabled", true);
		run("non-grav on restored");
		verifyReturn(baseline);
		QVERIFY(!QTest::currentTestFailed());
		settings.nonGrav()->setProperty("enabled", false);
		run("non-grav off restored");
		verifyReturn(gravity);
		QVERIFY(!QTest::currentTestFailed());
		settings.nonGrav()->setProperty("enabled", true);
		settings.weighting()->setProperty("debiasEnabled", false);
		run("debias off");
		QVERIFY(!QTest::currentTestFailed());
		settings.weighting()->setProperty("debiasEnabled", true);
		run("debias on restored");
		verifyReturn(baseline);
		QVERIFY(!QTest::currentTestFailed());
		if (manualSelection) {
			rows->resetManualSelection();
			run("manual selection reset");
			verifyReturn(automaticBaseline);
		}
	}

	void liveRecalculationUsesSnapshotAndPreservesDraftOnFailure()
	{
		if (!qEnvironmentVariableIsSet("ORBITAL_UI_LIVE_FIT"))
			QSKIP("Live fit is opt-in");
		const auto oldDirectory = getExecutableDir();
		auto* oldSource = getDefaultEphemerisSource();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source source;
		setDefaultEphemerisSource(&source);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(oldSource);
			getExecutableDir() = oldDirectory;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		settings.perturbations()->setProperty("enabled", true);
		settings.perturbations()->setUseSb441(false);
		settings.nonGrav()->setProperty("enabled", false);
		settings.weighting()->setProperty("useWeights", true);
		settings.weighting()->setProperty("debiasEnabled", true);
		settings.sigmaClip()->setProperty("enabled", true);
		QTemporaryDir inputs;
		const auto file = inputs.filePath("observations.xml");
		QVERIFY(QFile::copy(QStringLiteral(ORBITAL_UI_SOURCE_DIR) +
		                        "/../OrbitalCore/tests/data/2020VV_ADES_2022.ESA11_table_replay.xml",
		                    file));
		QVERIFY(session.openInput(QUrl::fromLocalFile(file)));
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
		QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
		QVERIFY(session.solution()->available());
		QVERIFY(QFile::remove(file));
		auto* rows = session.solution()->observations();
		const auto oldOrbit = session.solution()->snapshot()->fit.orbit;
		const auto oldTotal = rows->totalCount();
		const auto oldId = session.solution()->snapshot()->observations.all[1].id;
		bool sawIod = false;
		const auto observe = connect(&session, &ApplicationSession::changed, this, [&] {
			sawIod = sawIod || session.progress().contains("IOD") || session.progress().contains("ATT-EXT");
		});
		for (int i = 2; i < oldTotal; ++i)
			rows->setIncluded(i, false);
		const auto pending = rows->pendingCount();
		QVERIFY(pending > 0);
		session.run();
		QVERIFY(!rows->editingEnabled());
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 30000);
		QVERIFY(!session.error().isEmpty());
		QCOMPARE(rows->pendingCount(), pending);
		QCOMPARE(session.solution()->snapshot()->fit.orbit, oldOrbit);
		QVERIFY(rows->editingEnabled());
		rows->discardEdits();
		rows->setIncluded(1, false);
		settings.nonGrav()->setProperty("model", 4);
		settings.nonGrav()->setProperty("solveMask", 2);
		settings.nonGrav()->setProperty("enabled", true);
		session.run();
		settings.nonGrav()->setProperty("enabled", false);
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
		QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
		QVERIFY(session.solution()->nonGravEnabled());
		QCOMPARE(rows->pendingCount(), 0);
		QCOMPARE(rows->manualCount(), 1);
		QCOMPARE(rows->totalCount(), oldTotal);
		QCOMPARE(rows->overrides().front().id, oldId);
		QVERIFY(std::ranges::find(session.solution()->snapshot()->fit.keptObsIndices, 1u) ==
		        session.solution()->snapshot()->fit.keptObsIndices.end());
		rows->setIncluded(2, false);
		const auto beforeCancel = rows->pendingCount();
		session.run();
		session.cancel();
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 30000);
		QCOMPARE(session.progress(), QString("Cancelled"));
		QCOMPARE(rows->pendingCount(), beforeCancel);
		QVERIFY(session.solution()->nonGravEnabled());
		session.run();
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
		QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
		QVERIFY(!session.solution()->nonGravEnabled());
		QCOMPARE(session.solution()->snapshot()->fit.uncertainty.fitParameterCovariance.rows(), 6);
		QVERIFY(sawIod);
		disconnect(observe);
	}

	void equinoctialValuesAndCovarianceEpoch()
	{
		namespace equ = orbital_core::hazard::equinoctial;
		namespace cot = orbital_core::cometary;
		auto result = fixture();
		const auto& orbit = *result.fit.orbit;
		const auto r = orbit.getR(), v = orbit.getV() * Constants::GAUSS_K;
		equ::Vector6 state;
		state << r.x, r.y, r.z, v.x, v.y, v.z;
		const auto elements = equ::fromCartesian(state).value();
		equ::Matrix6 toState;
		QVERIFY(equ::toCartesian(elements, &toState));
		equ::Matrix6 expected = equ::Matrix6::Identity() * 1e-12;
		expected(0, 1) = expected(1, 0) = 3e-13;
		expected(1, 5) = expected(5, 1) = -2e-13;
		auto& uncertainty = result.fit.uncertainty;
		uncertainty.cartesianStateCovariance = (toState * expected * toState.transpose()).eval();
		uncertainty.cartesianStateCovarianceFormal = uncertainty.cartesianStateCovariance;
		SolutionViewModel model;
		auto verifyRows = [&] {
			const auto rows = model.equinoctialElements();
			QCOMPARE(rows.size(), 6);
			const QStringList names{"a", "h", "k", "p", "q", QString::fromUtf8("\xce\xbb")};
			for (int i = 0; i < 6; ++i) {
				const auto row = rows[i].toMap();
				const double scale = i == 5 ? Constants::RAD_TO_DEG : 1.0;
				QCOMPARE(row["label"].toString(), names[i]);
				QVERIFY(std::abs(row["value"].toDouble() - elements(i) * scale) < 1e-10);
				QVERIFY(std::abs(row["sigma"].toDouble() / (1e-6 * scale) - 1.0) < 1e-8);
			}
		};
		QVERIFY(model.publish(result, "equ"));
		verifyRows();

		const auto conventional = cot::fromCartesian(state).value();
		const auto cotJacobian = cot::cartesianToElementsJacobian(conventional).value();
		const equ::Matrix6 cotCovariance = cotJacobian * uncertainty.cartesianStateCovariance * cotJacobian.transpose();
		equ::Matrix6 toKeplerian = equ::Matrix6::Identity();
		toKeplerian(0, 0) = 1.0 / (1.0 - orbit.getExcentricity());
		toKeplerian(0, 1) = orbit.getPerigee() / std::pow(1.0 - orbit.getExcentricity(), 2);
		uncertainty.elementCovariance = (toKeplerian * cotCovariance * toKeplerian.transpose()).eval();
		uncertainty.hasElementCovariance = true;
		result.fit.fitEpochOrbit->setEpoch(orbit.getEpoch() - 10);
		uncertainty.cartesianStateCovariance *= 1e6;
		QVERIFY(model.publish(result, "equ"));
		verifyRows();

		uncertainty.hasElementCovariance = false;
		uncertainty.cometaryElementCovariance.emplace();
		uncertainty.cometaryElementCovariance->covariance = cotCovariance;
		uncertainty.cometaryElementCovariance->epochJdTdb = orbit.getEpoch();
		QVERIFY(model.publish(result, "equ"));
		verifyRows();
		uncertainty.cometaryElementCovariance->epochJdTdb -= 1;
		QVERIFY(model.publish(result, "equ"));
		for (const auto& row : model.equinoctialElements())
			QVERIFY(std::isnan(row.toMap()["sigma"].toDouble()));

		result.fit.orbit = make_orbit_from_states_vectors(
		    StateVectors(Vec3D(1, 0, 0), Vec3D(0, 0.04, 0) / Constants::GAUSS_K), orbit.getEpoch());
		QVERIFY(model.publish(result, "hyperbolic"));
		QVERIFY(model.equinoctialElements().isEmpty());
	}

	void observationSelectionFilteringAndSigns()
	{
		ObservationTableModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		model.setResult(fixture());
		QCOMPARE(model.totalCount(), 3);
		QCOMPARE(model.columnCount(), 11);
		QVERIFY(model.data(model.index(0, 4)).toString().startsWith("-0"));
		model.selectRow(0);
		model.sortColumn(1, true);
		QCOMPARE(model.selectedSourceRow(), 0);
		QCOMPARE(model.selectedProxyRow(), 2);
		QCOMPARE(model.sortedColumn(), 1);
		QVERIFY(model.sortDescending());
		QCOMPARE(model.data(model.index(0, 0)).toInt(), 3);
		model.setExcludedOnly(true);
		QCOMPARE(model.selectedProxyRow(), -1);
		QCOMPARE(model.rowCount(), 1);
		QCOMPARE(model.data(model.index(0, 0)).toInt(), 2);
		model.selectRow(0);
		QVERIFY(model.selectionDetails().contains("B02"));
		model.setFilterText("A01");
		QCOMPARE(model.rowCount(), 0);
		model.setExcludedOnly(false);
		QCOMPARE(model.rowCount(), 2);
		model.setResult(fixture());
		QCOMPARE(model.selectedSourceRow(), -1);
	}

	void equinoctialModeAvailability()
	{
		SolutionViewModel model;
		auto result = fixture();
		QVERIFY(model.publish(result, "fixture"));
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData("import QtQuick\nimport TrajectaUI\nOrbitElementsTable { width: 350 }", QUrl());
		std::unique_ptr<QObject> root(
		    component.createWithInitialProperties({{"solution", QVariant::fromValue(&model)}}));
		QVERIFY2(root, qPrintable(component.errorString()));
		auto* item = qobject_cast<QQuickItem*>(root.get());
		auto* mode = findItem(item, "orbit-elements-mode");
		auto* equButton = findItem(item, "orbit-elements-mode-option-1");
		QVERIFY(mode && equButton);
		QVERIFY(equButton->isEnabled());
		QVERIFY(QMetaObject::invokeMethod(mode, "activated", Q_ARG(int, 1)));
		QCOMPARE(mode->property("currentIndex").toInt(), 1);
		result.fit.orbit = make_orbit_from_states_vectors(
		    StateVectors(Vec3D(1, 0, 0), Vec3D(0, 0.04, 0) / Constants::GAUSS_K), 2460000.5);
		QVERIFY(model.publish(result, "hyperbolic"));
		QVERIFY(!equButton->isEnabled());
		QCOMPARE(mode->property("currentIndex").toInt(), 0);
		QVERIFY(QMetaObject::invokeMethod(mode, "activated", Q_ARG(int, 1)));
		QCOMPARE(mode->property("currentIndex").toInt(), 0);
		QVERIFY(model.publish(fixture(), "fixture"));
		QVERIFY(equButton->isEnabled());
	}

	void radarReferenceFrequencyDisplaysPhysicalRangeRate()
	{
		auto result = fixture();
		result.observations.all[0].radar =
		    orbital_core::obs::RadarObservation{.dopplerHz = orbital_core::obs::RadarMeasurement{-2e-8, 1e-10},
		                                        .frequencyHz = 1,
		                                        .frequencyIsReference = true,
		                                        .transmitter = "251"};
		ObservationTableModel model;
		model.setResult(result);
		model.selectSourceRow(0);
		QCOMPARE(model.headerData(4, Qt::Horizontal).toString(), QString("Range rate, m/s"));
		QCOMPARE(model.selection()["frequency"].toString(), QString("-"));
		const auto channels = model.selection()["radarChannels"].toList();
		QCOMPARE(channels.size(), 1);
		QCOMPARE(channels.front().toMap()["label"].toString(), QString("Range rate, m/s"));
		QVERIFY(std::abs(model.index(0, 4).data(Qt::UserRole + 1).toDouble() - 2.99792458) < 1e-9);
		QVERIFY(channels.front().toMap()["sigma"].toDouble() > 0);
	}

	void radarRowsHaveNativeUnitsAndStableManualSelection()
	{
		auto result = fixture();
		auto& observation = result.observations.all[0];
		observation.radar =
		    orbital_core::obs::RadarObservation{.delaySeconds = orbital_core::obs::RadarMeasurement{192.02850713, 4e-6},
		                                        .dopplerHz = orbital_core::obs::RadarMeasurement{-102512.91, 0.25},
		                                        .frequencyHz = 2380e6,
		                                        .transmitter = "251"};
		observation.RA = observation.DEC = std::numeric_limits<double>::quiet_NaN();
		auto residuals = result.fit.allResiduals.data();
		residuals[0].radar =
		    Residual::Radar{.delaySeconds = Residual::RadarChannel{192.02850713, 192.02850913, -2e-6, 4e-6},
		                    .dopplerHz = Residual::RadarChannel{-102512.91, -102512.81, -0.1, 0.25},
		                    .transmitter = "251",
		                    .frequencyHz = 2380e6,
		                    .chiSquared = 0.41};
		residuals[0].dRA_arcsec = residuals[0].dDEC_arcsec = std::numeric_limits<double>::quiet_NaN();
		result.fit.allResiduals.clear();
		for (const auto& residual : residuals)
			result.fit.allResiduals.add(residual);
		ObservationTableModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		model.setResult(result);
		QCOMPARE(model.totalCount(), 3);
		QCOMPARE(model.radarCount(), 1);
		QCOMPARE(model.rowCount(), 2);
		model.selectSourceRow(0);
		QVERIFY(model.radarMode());
		QCOMPARE(model.rowCount(), 1);
		QCOMPARE(model.selection()["ra"].toString(), QString("-"));
		QCOMPARE(model.selection()["frequency"].toString(), QString("2380.000000"));
		const auto channels = model.selection()["radarChannels"].toList();
		QCOMPARE(channels.size(), 2);
		QCOMPARE(channels[0].toMap()["residual"].toString(), QString("-2.000"));
		model.setIncluded(0, false);
		QCOMPARE(model.overrides().size(), 1u);
		QCOMPARE(model.overrides()[0].id, observation.id);
		model.setRadarMode(false);
		QCOMPARE(model.rowCount(), 2);
		QCOMPARE(model.pendingCount(), 1);
		model.setRadarMode(true);
		model.setResult(fixture());
		QVERIFY(!model.radarMode());
		QCOMPARE(model.rowCount(), 3);

		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QVERIFY(session.solution()->publish(result, "radar-test"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QVERIFY(QTest::qWaitForWindowExposed(window));
		window->setProperty("currentPage", 1);
		auto* observations = session.solution()->observations();
		observations->setRadarMode(true);
		observations->selectSourceRow(0);
		auto* page = findItem(window->contentItem(), "observations-page");
		QVERIFY(page);
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		for (const auto& theme : {QString("graphite"), QString("light")}) {
			session.setTheme(theme);
			window->resize(1512, 982);
			QTest::qWait(120);
			QVERIFY(QMetaObject::invokeMethod(page, "openDetails"));
			QTRY_VERIFY(findItem(window->contentItem(), "observations-inspector")->isVisible());
			QTest::qWait(120);
			if (!captures.isEmpty())
				QVERIFY(window->grabWindow().save(captures + "/radar-" + theme + ".png"));
			window->resize(1024, 720);
			QTest::qWait(120);
			QVERIFY(QMetaObject::invokeMethod(page, "openDetails"));
			auto* dialog = window->findChild<QObject*>("observations-page-details");
			QVERIFY(dialog);
			QTRY_VERIFY(dialog->property("visible").toBool());
			QVERIFY(dialog->property("height").toDouble() <= window->height() - 40);
			QTest::qWait(120);
			if (!captures.isEmpty())
				QVERIFY(window->grabWindow().save(captures + "/radar-dialog-" + theme + ".png"));
			QTest::keyClick(window, Qt::Key_Escape);
			QTRY_VERIFY(!dialog->property("visible").toBool());
		}
		window->close();
	}

	void observationMetadataAndWeights()
	{
		auto result = fixture();
		auto& observation = result.observations.all[0];
		observation.note2 = 'B';
		observation.catalogCode = 'X';
		observation.programCode = 'a';
		observation.adesProgramCode = "123";
		observation.biasRaArcsec = 0.123;
		observation.biasDecArcsec = -0.234;
		result.observations.preparation = orbital_core::ObservationPreparation{.debiasRequested = true};
		auto residuals = result.fit.allResiduals.data();
		residuals[0].astrometricWeight = ResidualWeight{.sigmaRaArcsec = 0.25,
		                                                .sigmaDecArcsec = 0.5,
		                                                .chiSquared = 25.0,
		                                                .weighted = true,
		                                                .postFitChiSquared = 36.0};
		result.fit.allResiduals.clear();
		for (const auto& residual : residuals)
			result.fit.allResiduals.add(residual);
		ObservationTableModel model;
		model.setResult(result);
		model.selectSourceRow(0);
		const auto fields = model.selection();
		QCOMPARE(fields["technique"].toString(), QString("CMOS"));
		QCOMPARE(fields["catalog"].toString(), QString("Gaia EDR3 (X)"));
		QCOMPARE(fields["program"].toString(), QString("123"));
		QCOMPARE(fields["sigmaRa"].toString(), QString("0.250"));
		QCOMPARE(fields["sigmaDec"].toString(), QString("0.500"));
		QCOMPARE(fields["biasRa"].toString(), QString("0.123"));
		QCOMPARE(fields["biasDec"].toString(), QString("-0.234"));
		QCOMPARE(fields["chi"].toString(), QString("5.00"));
		QCOMPARE(fields["postFitChi"].toString(), QString("6.00"));
		QCOMPARE(model.data(model.index(0, 7)).toString(), QString("0.250"));
		QCOMPARE(model.data(model.index(0, 8)).toString(), QString("0.500"));
		QCOMPARE(model.data(model.index(0, 9)).toString(), QString("5.00"));
		QCOMPARE(model.data(model.index(0, 10)).toString(), QString("Included"));
		QVERIFY(!model.data(model.index(1, 7)).isValid());
		model.setFilterText("123");
		QCOMPARE(model.count(), 1);
		model.setFilterText("gaia edr3");
		QCOMPARE(model.count(), 1);
		model.setFilterText("CMOS");
		QCOMPARE(model.count(), 1);
		model.setStationFilter("B02");
		model.setExcludedOnly(true);
		model.clearFilters();
		QCOMPARE(model.count(), 3);
		QVERIFY(model.filterText().isEmpty());
		QVERIFY(model.stationFilter().isEmpty());
		QVERIFY(!model.excludedOnly());
		QCOMPARE(model.selectedSourceRow(), 0);
		QVERIFY(fields["debiasEnabled"].toBool());
		residuals[0].astrometricWeight->weighted = false;
		result.fit.allResiduals.clear();
		for (const auto& residual : residuals)
			result.fit.allResiduals.add(residual);
		model.setResult(result);
		model.selectSourceRow(0);
		QCOMPARE(model.selection()["weighting"].toString(), QString("Unweighted"));
		QCOMPARE(model.selection()["sigmaRa"].toString(), QString("-"));
		QCOMPARE(model.selection()["chi"].toString(), QString("-"));
		QVERIFY(!model.data(model.index(0, 7)).isValid());
		QVERIFY(!model.data(model.index(0, 9)).isValid());
		model.selectSourceRow(1);
		QCOMPARE(model.selection()["technique"].toString(), QString("Not supplied"));
		QCOMPARE(model.selection()["catalog"].toString(), QString("Not supplied"));
		QCOMPARE(model.selection()["sigmaDec"].toString(), QString("-"));
	}

	void cometClassesShownSeparately()
	{
		auto result = fixture();
		result.observations.classificationDesignation = "2P";
		auto& orbit = *result.fit.orbit;
		orbit.setType(OrbitType::Elliptical);
		orbit.setSemiMajorAxis(2.22);
		orbit.setExcentricity(0.847);
		orbit.setInclination(Convertors::degreesToRadians(11.3));
		orbit.setOrbitalPeriod(1210.0);
		SolutionViewModel solution;
		QVERIFY(solution.publish(result, "test"));
		QCOMPARE(solution.orbit()->orbitalClasses(),
		         QStringList({"Comet", "Jupiter family (period)", "Encke type (Tj)"}));
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(
		    &engine, QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/SolutionHeader.qml"));
		QVERIFY2(component.isReady(), qPrintable(component.errorString()));
		std::unique_ptr<QObject> header(component.createWithInitialProperties(
		    {{"solution", QVariant::fromValue(&solution)}, {"width", 1024}, {"height", 73}, {"viewportWidth", 1024}}));
		QVERIFY2(header != nullptr, qPrintable(component.errorString()));
		auto* headerItem = qobject_cast<QQuickItem*>(header.get());
		QVERIFY(headerItem);
		QQuickWindow window;
		window.resize(1024, 200);
		headerItem->setParentItem(window.contentItem());
		window.show();
		const auto countBadges = [](QQuickItem* root, const auto& self) -> int {
			int count = root->objectName() == "solution-orbital-class" ? 1 : 0;
			for (auto* child : root->childItems())
				count += self(child, self);
			return count;
		};
		QTRY_COMPARE(countBadges(headerItem, countBadges), 3);
		result.fit.converged = false;
		QVERIFY(solution.publish(result, "test"));
		QCOMPARE(solution.orbit()->orbitalClasses(), QStringList({"Comet"}));
		QTRY_COMPARE(countBadges(headerItem, countBadges), 1);
		headerItem->setParentItem(nullptr);
		window.close();
	}

	void observationTimeAndReplacedRecordPresentation()
	{
		QCOMPARE(FormatUtils::dec(Convertors::degreesToRadians(22.0 + 1.0 / 60.0)), QString("+22 01 00.0"));
		QCOMPARE(FormatUtils::dec(Convertors::degreesToRadians(-(22.0 + 1.0 / 60.0))), QString("-22 01 00.0"));
		QCOMPARE(FormatUtils::dec(-0.0), QString("-0 00 00.0"));
		QCOMPARE(FormatUtils::ra((24.0 - 0.001 / 3600.0) * Constants::PI / 12.0), QString(" 0:00:00.00"));

		auto result = fixture();
		auto& observation = result.observations.all[0];
		observation.note2 = 'X';
		observation.timeScale = orbital_core::time::Scale::UT1;
		observation.DEC = Convertors::degreesToRadians(22.0 + 1.0 / 60.0);
		result.fit.keptObsIndices = {2};
		auto residuals = result.fit.allResiduals.data();
		residuals[0].astrometricWeight = ResidualWeight{.sigmaRaArcsec = 999.0,
		                                                .sigmaDecArcsec = 999.0,
		                                                .chiSquared = 0.0036,
		                                                .weighted = true,
		                                                .postFitChiSquared = 0.0036};
		result.fit.allResiduals.clear();
		for (const auto& residual : residuals)
			result.fit.allResiduals.add(residual);

		ObservationTableModel model;
		model.setResult(result);
		model.selectSourceRow(0);
		QCOMPARE(model.headerData(1, Qt::Horizontal).toString(), QString("Date and time"));
		QVERIFY(model.data(model.index(0, 1)).toString().endsWith(" UT1"));
		QCOMPARE(model.data(model.index(0, 4)).toString(), QString("+22 01 00.0"));
		QVERIFY(!model.data(model.index(0, 7)).isValid());
		QVERIFY(!model.data(model.index(0, 8)).isValid());
		QVERIFY(!model.data(model.index(0, 9)).isValid());
		QCOMPARE(model.data(model.index(0, 10)).toString(), QString("Replaced"));
		QCOMPARE(model.selection()["timeScale"].toString(), QString("UT1"));
		QCOMPARE(model.selection()["weighting"].toString(), QString("Replaced by MPC"));
		QCOMPARE(model.selection()["sigmaRa"].toString(), QString("-"));
		QCOMPARE(model.selection()["chi"].toString(), QString("-"));
		QVERIFY(model.selectionDetails().contains(" UT1"));
		QVERIFY(model.data(model.index(1, 1)).toString().endsWith(" UTC"));

		result.fit.keptObsIndices = {0, 2};
		model.setResult(result);
		model.selectSourceRow(0);
		QCOMPARE(model.selection()["weighting"].toString(), QString("Weighted"));
		QCOMPARE(model.selection()["sigmaRa"].toString(), QString("999.000"));
		QCOMPARE(model.data(model.index(0, 10)).toString(), QString("Included"));
	}

	void observationPageInteractions()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto result = fixture();
		const auto seed = result.observations.all.front();
		result.observations.all.clear();
		result.observations.format = "ADES XML";
		result.fit.keptObsIndices.clear();
		result.fit.allResiduals.clear();
		result.fit.residuals.clear();
		for (int i = 0; i < 120; ++i) {
			auto row = seed;
			row.id = orbital_core::newObservationId();
			row.JD += i * 0.01;
			row.obsCode = i % 2 ? "B02" : "A01";
			row.note2 = 'C';
			row.catalogCode = 'V';
			row.astCat = "Gaia2";
			row.adesProgramCode = i % 2 ? "P12" : "P34";
			result.observations.all.push_back(row);
			Residual residual{};
			residual.dRA_arcsec = (i % 7 - 3) * 0.123;
			residual.dDEC_arcsec = (i % 5 - 2) * 0.132;
			residual.astrometricWeight = ResidualWeight{.sigmaRaArcsec = 0.2 + 0.01 * (i % 10),
			                                            .sigmaDecArcsec = 0.3,
			                                            .chiSquared = 0.5 * (i % 10),
			                                            .weighted = true,
			                                            .postFitChiSquared = 0.6 * (i % 10)};
			result.fit.allResiduals.add(residual);
			if (i % 5) {
				result.fit.keptObsIndices.push_back(i);
				result.fit.residuals.add(residual);
			}
		}
		QVERIFY(session.solution()->publish(result, "test"));
		auto* observations = session.solution()->observations();
		QCOMPARE(observations->sourceFormat(), QString("ADES XML"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QVERIFY(QTest::qWaitForWindowExposed(window));
		window->resize(1512, 982);
		window->setProperty("currentPage", 1);
		const auto click = [&](const QString& name) {
			auto* item = findItem(window->contentItem(), name);
			QVERIFY2(item && item->isVisible(), qPrintable(name));
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
			                  item->mapToScene({item->width() / 2, item->height() / 2}).toPoint());
		};
		auto* page = findItem(window->contentItem(), "observations-page");
		auto* inspector = findItem(window->contentItem(), "observations-inspector");
		auto* table = findItem(window->contentItem(), "observations-table");
		auto* view = findItem(window->contentItem(), "observations-data-view");
		QVERIFY(page && inspector && table && view);
		QTRY_VERIFY(page->property("inlineInspector").toBool());
		QVERIFY(!inspector->isVisible());
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			observations->setIncluded(1, false);
			for (const QSize size : {QSize(1024, 720), QSize(1512, 982), QSize(2554, 1312)}) {
				window->resize(size);
				QTest::qWait(120);
				const auto bounds = table->mapRectToScene(table->boundingRect());
				QVERIFY(bounds.width() > 850);
				QVERIFY(bounds.height() > 300);
				QVERIFY(bounds.right() <= window->width());
				QVERIFY(bounds.bottom() <= window->height());
				if (!captures.isEmpty())
					QVERIFY(window->grabWindow().save(captures + "/observations-" + theme + "-" +
					                                  QString::number(size.width()) + ".png"));
			}
		}
		observations->discardEdits();
		window->resize(1512, 982);
		session.setTheme("graphite");
		QTest::qWait(120);
		click("observation-included-0");
		QCOMPARE(observations->pendingCount(), 1);
		QVERIFY(!inspector->isVisible());
		QCOMPARE(observations->overrides().front().id, result.observations.all[0].id);
		click("observations-discard");
		QCOMPARE(observations->pendingCount(), 0);
		observations->sortColumn(1, true);
		observations->setStationFilter("B02");
		QTest::qWait(100);
		click("observation-included-119");
		QCOMPARE(observations->pendingCount(), 1);
		QCOMPARE(observations->overrides().front().id, result.observations.all[119].id);
		observations->setEditingEnabled(false);
		QVERIFY(!findItem(window->contentItem(), "observation-included-119")->isEnabled());
		observations->setEditingEnabled(true);
		click("observations-discard");
		observations->clearFilters();
		observations->sortColumn(0, false);
		QTest::qWait(100);
		click("observation-cell-0-1");
		QCOMPARE(observations->selectedSourceRow(), 0);
		QTRY_VERIFY(inspector->isVisible());
		QTest::qWait(100);
		QVERIFY(table->mapRectToScene(table->boundingRect()).right() < inspector->mapToScene({}).x());
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/observations-inspector.png"));
		click("observations-header-1");
		click("observations-header-1");
		QVERIFY(observations->sortDescending());
		QCOMPARE(observations->selectedSourceRow(), 0);
		QCOMPARE(observations->selectedProxyRow(), 119);
		table->forceActiveFocus();
		QTest::keyClick(window, Qt::Key_Home);
		QCOMPARE(observations->selectedSourceRow(), 119);
		QTest::keyClick(window, Qt::Key_Down);
		QCOMPARE(observations->selectedSourceRow(), 118);
		QTest::keyClick(window, Qt::Key_End);
		QCOMPARE(observations->selectedSourceRow(), 0);
		click("observations-inspector-close");
		QTRY_VERIFY(!inspector->isVisible());
		observations->setFilterText("no matches");
		QCOMPARE(observations->count(), 0);
		QTRY_VERIFY(findItem(page, "observations-empty")->isVisible());
		click("observations-empty-action");
		QCOMPARE(observations->count(), 120);
		click("observations-excluded-filter");
		QCOMPARE(observations->count(), 24);
		observations->setFilterText("P12");
		QCOMPARE(observations->count(), 12);
		click("observations-clear-search");
		QCOMPARE(observations->count(), 24);
		observations->clearFilters();
		observations->sortColumn(0, false);
		window->resize(1024, 720);
		QTest::qWait(150);
		click("observation-cell-1-1");
		auto* dialog = window->findChild<QObject*>("observations-page-details");
		QVERIFY(dialog);
		QTRY_VERIFY(dialog->property("visible").toBool());
		QCOMPARE(observations->selectedSourceRow(), 1);
		QVERIFY(dialog->property("height").toDouble() <= window->height() - 40);
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/observations-compact-details.png"));
		QTest::keyClick(window, Qt::Key_Escape);
		QTRY_VERIFY(!dialog->property("visible").toBool());
		window->close();
	}

	void observationCatalogs()
	{
		auto result = fixture();
		auto& observation = result.observations.all[0];
		ObservationTableModel model;
		const auto check = [&](char code, const std::string& astCat, const QString& expected) {
			observation.catalogCode = code;
			observation.astCat = astCat;
			model.setResult(result);
			model.selectSourceRow(0);
			QCOMPARE(model.selection()["catalog"].toString(), expected);
		};
		check('X', "GaIa3E", "Gaia EDR3 (X)");
		check('7', "ZZCAT", "ZZCAT (7)");
		check('&', "SDSS9", "SDSS DR9 (&)");
		check('9', "", "UBSC (9)");
		check('9', "UBSC", "UBSC (9)");
		check(' ', "FutureCatalog_DR9", "FutureCatalog_DR9 (unrecognized)");
		check(' ', "UNK", "Unknown");
		check(' ', "", "Not supplied");
		check('\0', "", "Not supplied");
		check('?', "", "MPC code ?");
	}

	void photometryAndMoidSnapshot()
	{
		auto result = fixture();
		result.fit.allResiduals.clear();
		Residual photometry{};
		photometry.mag = 20;
		photometry.magBand = 'V';
		photometry.helioDistAu = 1;
		photometry.geoDistAu = 1;
		result.fit.allResiduals.add(photometry);
		orbital_core::hazard::PlanetMoidTable moids;
		moids[orbital_core::hazard::Planet::Earth] = 0.001;
		SolutionViewModel model;
		QVERIFY(model.publish(result, "saved", moids));
		QCOMPARE(model.absoluteMagnitude(), 20.0);
		QCOMPARE(model.magnitudeCount(), 1);
		const auto sizes = model.sizeEstimates();
		QCOMPARE(sizes.size(), 5);
		const auto darkest = sizes[0].toMap();
		QCOMPARE(darkest["label"].toString(), QString::fromUtf8("0.05\xe2\x80\x93"
		                                                        "0.10"));
		QCOMPARE(darkest["value"].toString(), QString::fromUtf8("420\xe2\x80\x93"
		                                                        "594"));
		QCOMPARE(darkest["unit"].toString(), QString("m"));
		QCOMPARE(darkest["diameterHiKm"].toDouble(), 1329.0 / std::sqrt(0.05) * 1e-4);
		QCOMPARE(sizes[4].toMap()["value"].toString(), QString::fromUtf8("133\xe2\x80\x93"
		                                                                 "188"));
		QCOMPARE(model.typicalSize()["value"].toString(), QString::fromUtf8("266\xe2\x80\x93"
		                                                                    "594"));
		QCOMPARE(model.earthMoidKm(), 0.001 * (Constants::AU_IN_METERS / 1000.0));
		QCOMPARE(model.planetMoids().size(), 8);
		QCOMPARE(model.planetMoids()[2].toMap()["planet"].toString(), QString("Earth"));
		QVERIFY(std::isnan(model.planetMoids()[0].toMap()["au"].toDouble()));
		result.fit.cancelled = true;
		QVERIFY(!model.publish(result, "cancelled"));
		QCOMPARE(model.absoluteMagnitude(), 20.0);
		QCOMPARE(model.earthMoidKm(), 0.001 * (Constants::AU_IN_METERS / 1000.0));
		QVERIFY(model.publish(fixture(), "no-photometry"));
		QVERIFY(std::isnan(model.absoluteMagnitude()));
		QVERIFY(model.sizeEstimates().isEmpty());
		QVERIFY(model.typicalSize().isEmpty());
		QVERIFY(std::isnan(model.earthMoidKm()));
	}

	void nonGravResultSnapshot()
	{
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		QSettings().clear();
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		settings.perturbations()->setProperty("enabled", false);
		auto* model = session.solution();
		QVERIFY(!model->nonGravEnabled());
		QVERIFY(model->nonGravParameters().isEmpty());
		auto result = fixture();
		QVERIFY(model->publish(result, "off"));
		QCOMPARE(model->nonGravModel(), QString("Off"));
		for (const auto& row : model->nonGravParameters()) {
			QCOMPARE(row.toMap()["status"].toString(), QString("Not used"));
			QVERIFY(std::isnan(row.toMap()["value"].toDouble()));
		}
		auto& ng = result.fit.finalNonGrav;
		ng.setEnabled(true);
		ng.setModel(NgModel::Yarkovsky);
		ng.setParams(0, -1.234567e-13, 0);
		ng.setDT(12);
		auto& covariance = result.fit.uncertainty.fitParameterCovariance;
		covariance = orbital_core::FitParameterCovariance::Identity(7, 7) * 1e-8;
		covariance(6, 6) = 4e-28;
		result.fit.uncertainty.fitParameterCovarianceFormal = orbital_core::FitParameterCovariance::Identity(7, 7);
		QVERIFY(model->publish(result, "yarkovsky"));
		QVERIFY(model->nonGravEnabled());
		QCOMPARE(model->nonGravModel(), QString("Yarkovsky"));
		const auto rows = model->nonGravParameters();
		QCOMPARE(rows.size(), 4);
		QCOMPARE(rows[1].toMap()["value"].toDouble(), -1.234567e-13);
		QCOMPARE(rows[1].toMap()["sigma"].toDouble(), 2e-14);
		QCOMPARE(rows[1].toMap()["status"].toString(), QString("Estimated"));
		QCOMPARE(rows[1].toMap()["unit"].toString(), QStringLiteral("AU/day\u00b2"));
		for (int i : {0, 2, 3}) {
			QVERIFY(!rows[i].toMap()["used"].toBool());
			QVERIFY(std::isnan(rows[i].toMap()["value"].toDouble()));
			QVERIFY(std::isnan(rows[i].toMap()["sigma"].toDouble()));
		}
		auto* draft = session.calculationSettings();
		draft->reset();
		draft->nonGrav()->setProperty("model", 1);
		draft->nonGrav()->setProperty("solveMask", 15);
		QVERIFY(session.applySettings());
		QCOMPARE(model->nonGravModel(), QString("Yarkovsky"));
		QCOMPARE(model->nonGravParameters()[1].toMap()["value"].toDouble(), -1.234567e-13);
		ng.setEnabled(false);
		result.fit.cancelled = true;
		QVERIFY(!model->publish(result, "cancelled"));
		QVERIFY(!model->publish({}, "failed"));
		QCOMPARE(model->nonGravModel(), QString("Yarkovsky"));
		QCOMPARE(model->nonGravParameters()[1].toMap()["sigma"].toDouble(), 2e-14);
		QVERIFY(model->publish(fixture(), "off-again"));
		QVERIFY(!model->nonGravEnabled());
	}

	void nonGravMasksAndUncertainty()
	{
		SolutionViewModel model;
		auto result = fixture();
		auto& ng = result.fit.finalNonGrav;
		ng.setEnabled(true);
		ng.setParams(1.25e-8, -2.5e-9, 3.75e-10);
		ng.setDT(12.5);
		const double values[]{1.25e-8, -2.5e-9, 3.75e-10, 12.5};
		for (auto type : {NgModel::MarsdenH2O, NgModel::MarsdenCO, NgModel::InverseSquare}) {
			ng.setModel(type);
			const int count = ng.numFitParams();
			auto& covariance = result.fit.uncertainty.fitParameterCovariance;
			covariance = orbital_core::FitParameterCovariance::Identity(6 + count, 6 + count);
			for (int i = 0; i < count; ++i)
				covariance(6 + i, 6 + i) = (i + 1) * (i + 1) * 1e-20;
			QVERIFY(model.publish(result, "estimated"));
			for (int i = 0; i < count; ++i) {
				const auto row = model.nonGravParameters()[i].toMap();
				QCOMPARE(row["value"].toDouble(), values[i]);
				QCOMPARE(row["status"].toString(), QString("Estimated"));
				QVERIFY(std::abs(row["sigma"].toDouble() - (i + 1) * 1e-10) < 1e-24);
			}
			QCOMPARE(model.nonGravParameters()[3].toMap()["unit"].toString(), QString("days"));
			ng.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2);
			covariance = orbital_core::FitParameterCovariance::Identity(8, 8);
			QVERIFY(model.publish(result, "partial"));
			const auto rows = model.nonGravParameters();
			QCOMPARE(rows[2].toMap()["status"].toString(), QString("Fixed"));
			QCOMPARE(rows[2].toMap()["value"].toDouble(), values[2]);
			QVERIFY(std::isnan(rows[2].toMap()["sigma"].toDouble()));
			QCOMPARE(rows[3].toMap()["status"].toString(), QString("Fixed"));
			ng.setFitParameterMask(NgParameterMask::None);
			QVERIFY(model.publish(result, "fixed"));
			QCOMPARE(model.nonGravParameters()[0].toMap()["status"].toString(), QString("Fixed"));
		}
		ng.setModel(NgModel::Yarkovsky);
		auto& covariance = result.fit.uncertainty.fitParameterCovariance;
		ng.setFitParameterMask(NgParameterMask::None);
		QVERIFY(model.publish(result, "fixed-yarkovsky"));
		QCOMPARE(model.nonGravParameters()[0].toMap()["status"].toString(), QString("Fixed"));
		QCOMPARE(model.nonGravParameters()[0].toMap()["value"].toDouble(), values[0]);
		ng.setFitParameterMask(NgParameterMask::A2);
		for (int invalid = 0; invalid < 6; ++invalid) {
			covariance = orbital_core::FitParameterCovariance::Identity(7, 7);
			switch (invalid) {
				case 0:
					covariance.resize(0, 0);
					break;
				case 1:
					covariance.conservativeResize(6, 6);
					break;
				case 2:
					covariance(6, 6) = -1;
					break;
				case 3:
					covariance(6, 6) = std::numeric_limits<double>::infinity();
					break;
				case 4:
					covariance(0, 6) = std::numeric_limits<double>::quiet_NaN();
					break;
				case 5:
					covariance.setZero();
					break;
			}
			QVERIFY(model.publish(result, "unavailable-sigma"));
			const auto row = model.nonGravParameters()[1].toMap();
			QCOMPARE(row["status"].toString(), QString("Estimated"));
			QCOMPARE(row["value"].toDouble(), values[1]);
			QVERIFY(std::isnan(row["sigma"].toDouble()));
		}
	}

	void nonGravPresentation()
	{
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		QSettings().clear();
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		session.setTheme("observatory");
		auto result = fixture();
		auto& ng = result.fit.finalNonGrav;
		ng.setEnabled(true);
		ng.setModel(NgModel::Yarkovsky);
		ng.setParams(0, -1.234567e-13, 0);
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(7, 7) * 4e-28;
		QVERIFY(session.solution()->publish(result, "fixture"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QTRY_VERIFY(window->isVisible());
		auto* summary = findItem(window->contentItem(), "nominal-non-grav");
		auto* dialog = window->findChild<QObject*>("orbit-details-dialog");
		QVERIFY(summary && dialog);
		QVERIFY(!findItem(summary, "non-grav-row-A1"));
		auto* a2 = findItem(summary, "non-grav-row-A2");
		QVERIFY(a2);
		QCOMPARE(findItem(a2, "non-grav-value")->property("text").toString(), QString("-1.234567e-13"));
		QCOMPARE(findItem(a2, "non-grav-sigma")->property("text").toString(), QString("2.000e-14"));
		const auto captureDir = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captureDir.isEmpty())
			QVERIFY(QDir().mkpath(captureDir));
		for (const QSize size : {QSize(1512, 982), QSize(1024, 720)}) {
			window->resize(size);
			QTest::qWait(150);
			auto* value = findItem(a2, "non-grav-value");
			QVERIFY(value->width() >= value->property("contentWidth").toDouble());
			if (!captureDir.isEmpty())
				QVERIFY(window->grabWindow().save(captureDir + "/non-grav-solution-" + QString::number(size.width()) +
				                                  ".png"));
			QVERIFY(QMetaObject::invokeMethod(dialog, "open"));
			QTest::qWait(150);
			auto* content = dialog->property("contentItem").value<QQuickItem*>();
			auto* details = findItem(content, "detailed-non-grav");
			QVERIFY(details);
			auto* sizeSection = findItem(content, "estimated-size");
			QVERIFY(sizeSection && !sizeSection->isVisible());
			auto* derived = findItem(content, "derived-row-q");
			auto* state = findItem(content, "state-row-X");
			QVERIFY(derived && state);
			QCOMPARE(derived->height(), state->height());
			QCOMPARE(derived->mapToItem(content, QPointF()).y(), state->mapToItem(content, QPointF()).y());
			for (const auto& field : {std::pair{derived, "derived-value"}, std::pair{state, "state-value"}}) {
				auto* value = findItem(field.first, field.second);
				QVERIFY(value);
				QVERIFY(value->width() >= value->property("contentWidth").toDouble());
			}
			QVERIFY(dialog->property("height").toDouble() <= window->height() - 56);
			auto* viewport = content->property("contentItem").value<QQuickItem*>();
			QVERIFY(viewport);
			auto* scrollbar = findItem(content, "orbit-details-scrollbar");
			QVERIFY(scrollbar);
			QCOMPARE(scrollbar->height(), content->property("availableHeight").toDouble());
			QCOMPARE(scrollbar->x() + scrollbar->width(), content->width());
			viewport->setProperty("contentY",
			                      std::max(0.0, viewport->property("contentHeight").toDouble() - viewport->height()));
			QTest::qWait(50);
			for (const QString label : {"A1", "A2", "A3", "DT"}) {
				auto* row = findItem(details, "non-grav-row-" + label);
				QVERIFY(row);
				QCOMPARE(row->height(), derived->height());
				auto* value = findItem(row, "non-grav-value");
				auto* sigma = findItem(row, "non-grav-sigma");
				QVERIFY(value && sigma);
				QVERIFY(value->width() >= value->property("contentWidth").toDouble());
				QVERIFY(sigma->width() >= sigma->property("contentWidth").toDouble());
				QVERIFY(value->x() + value->width() <= sigma->x());
				QVERIFY(value->mapToItem(row, {value->width(), 0}).x() < row->width() * 0.55);
				auto* status = findItem(row, "non-grav-status");
				QVERIFY(status->isVisible());
				QVERIFY(status->width() >= status->property("contentWidth").toDouble());
				QCOMPARE(status->property("text").toString(),
				         label == "A2" ? QString("Estimated") : QString("Not used"));
				QVERIFY(row->mapToItem(viewport, QPointF(0, row->height())).y() <= viewport->height() + 1);
			}
			if (!captureDir.isEmpty())
				QVERIFY(window->grabWindow().save(captureDir + "/non-grav-dialog-" + QString::number(size.width()) +
				                                  ".png"));
			QVERIFY(QMetaObject::invokeMethod(dialog, "close"));
		}
		ng.setModel(NgModel::MarsdenH2O);
		ng.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2);
		ng.setParams(1.234567e-8, -2.345678e-9, 3.456789e-10);
		ng.setDT(12.5);
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(8, 8) * 1e-20;
		QVERIFY(session.solution()->publish(result, "partial"));
		QTest::qWait(150);
		auto* a3 = findItem(summary, "non-grav-row-A3");
		QVERIFY(a3);
		QCOMPARE(findItem(a3, "non-grav-sigma")->property("text").toString(), QString("Fixed"));
		QVERIFY(QMetaObject::invokeMethod(dialog, "open"));
		QTest::qWait(150);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/non-grav-dialog-fixed.png"));
		QVERIFY(QMetaObject::invokeMethod(dialog, "close"));
		ng.setEnabled(false);
		QVERIFY(session.solution()->publish(result, "off"));
		QTest::qWait(150);
		QCOMPARE(findItem(summary, "non-grav-model")->property("text").toString(), QString("Off"));
		QVERIFY(!findItem(summary, "non-grav-row-A2"));
		window->close();
	}

	void covarianceUnitsAndCorrelation()
	{
		CovarianceTableModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		auto result = fixture();
		result.fit.uncertainty.fitParameterCovariance(0, 1) = -2e-9;
		result.fit.uncertainty.fitParameterCovariance(1, 0) = -2e-9;
		model.setResult(result.fit);
		QCOMPARE(model.dimension(), 6);
		QVERIFY(model.equinoctialAvailable());
		QVERIFY(model.equinoctial());
		QCOMPARE(model.epoch(), QString("60000.000000000 MJD TDB"));
		model.setEquinoctial(false);
		QCOMPARE(model.headerData(3, Qt::Horizontal, Qt::DisplayRole).toString(), QString("Vx"));
		QVERIFY(std::abs(model.data(model.index(0, 1), Qt::UserRole).toDouble() + 0.2) < 1e-12);
		model.setCorrelation(false);
		QCOMPARE(model.data(model.index(0, 1), Qt::UserRole).toDouble(), -2e-9);
		model.selectCell(0, 1);
		QCOMPARE(model.selection()["covariance"].toDouble(), -2e-9);
		QCOMPARE(model.selection()["rowUnit"].toString(), QString("AU"));
		QCOMPARE(model.selection()["rowSigma"].toDouble(), 1e-4);
		model.selectCell(-1, 100);
		QCOMPARE(model.selectedColumn(), 1);
		model.setEquinoctial(true);
		const auto equ = orbital_core::uncertainty::describeCovariance(
		    result.fit, orbital_core::uncertainty::CovarianceBasis::Equinoctial);
		QVERIFY(equ);
		QVERIFY(std::abs(model.data(model.index(5, 5), Qt::UserRole).toDouble() / equ->covariance(5, 5) -
		                 Constants::RAD_TO_DEG * Constants::RAD_TO_DEG) < 1e-10);
		QCOMPARE(model.weakDirection().size(), 6);
		QCOMPARE(model.weakSigma(), equ->spectrum->signedSigmas(5));
		result.fit.finalNonGrav.setModel(NgModel::Yarkovsky);
		result.fit.finalNonGrav.setEnabled(true);
		result.fit.uncertainty.fitParameterCovariance = orbital_core::FitParameterCovariance::Identity(7, 7);
		model.setResult(result.fit);
		QCOMPARE(model.dimension(), 7);
		QCOMPARE(model.headerData(6, Qt::Horizontal, Qt::DisplayRole).toString(), QString("A2"));
		model.selectCell(0, 6);
		QCOMPARE(model.selection()["columnUnit"].toString(), QStringLiteral("AU/day\u00b2"));
		model.setResult({});
		QCOMPARE(model.dimension(), 0);
		QVERIFY(model.epoch().isEmpty());
		QVERIFY(model.selection().isEmpty());
		QVERIFY(model.weakDirection().isEmpty());
	}

	void covarianceLayoutAndSelection()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QVERIFY(session.solution()->publish(fixture(), "test"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		auto* page = window->findChild<QObject*>("solution-page");
		QVERIFY(page);
		page->setProperty("currentTab", 1);
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		auto* model = session.solution()->covariance();
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			for (const QSize size : {QSize(1024, 720), QSize(1512, 982), QSize(2554, 1312)}) {
				window->resize(size);
				QTest::qWait(120);
				auto* panel = findItem(window->contentItem(), "covariance-panel");
				auto* matrix = findItem(window->contentItem(), "covariance-matrix");
				auto* inspector = findItem(window->contentItem(), "covariance-inspector");
				auto* cell = findItem(window->contentItem(), "covariance-cell-0-2");
				QVERIFY(panel && matrix && inspector && cell);
				QVERIFY(panel->x() >= 20);
				QVERIFY(panel->width() <= panel->parentItem()->width() - 40);
				QCOMPARE(panel->y(), 24.0);
				QCOMPARE(panel->property("stacked").toBool(), size.width() == 1024);
				QCOMPARE(matrix->property("contentWidth").toDouble(), matrix->width());
				QVERIFY(cell->width() >= 84);
				QCOMPARE(cell->height(), 68.0);
				const auto matrixBounds = matrix->mapRectToScene(matrix->boundingRect());
				const auto detailBounds = inspector->mapRectToScene(inspector->boundingRect());
				if (panel->property("stacked").toBool())
					QVERIFY(detailBounds.top() > matrixBounds.bottom());
				else
					QVERIFY(detailBounds.left() > matrixBounds.right());
				QVERIFY(detailBounds.right() <= panel->mapToScene({panel->width(), 0}).x() + 1);
				const auto click = [&](QQuickItem* item) {
					QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
					                  item->mapToScene({item->width() / 2, item->height() / 2}).toPoint());
				};
				click(cell);
				QCOMPARE(model->selectedRow(), 0);
				QCOMPARE(model->selectedColumn(), 2);
				QVERIFY(cell->property("checked").toBool());
				QCOMPARE(findItem(window->contentItem(), "covariance-pair-title")->property("text").toString(),
				         QString("a / k"));
				for (int mode : {0, 1}) {
					auto* button = findItem(window->contentItem(), "covariance-mode-option-" + QString::number(mode));
					QVERIFY(button);
					click(button);
					QCOMPARE(model->correlation(), mode == 0);
					const auto inspect = [&](const auto& self, QQuickItem* item) -> void {
						if (item->objectName().startsWith("covariance-cell-") && item->objectName().endsWith("-text"))
							QVERIFY(!item->property("truncated").toBool());
						for (auto* child : item->childItems())
							self(self, child);
					};
					inspect(inspect, matrix);
					if (!captures.isEmpty())
						QVERIFY(window->grabWindow().save(captures + "/covariance-" + theme + "-" +
						                                  QString::number(size.width()) + "-" + QString::number(mode) +
						                                  ".png"));
				}
			}
		}
		window->resize(1512, 982);
		auto full = fixture();
		full.fit.finalNonGrav.setModel(NgModel::MarsdenH2O);
		full.fit.finalNonGrav.setEnabled(true);
		const int count = 6 + full.fit.finalNonGrav.numFitParams();
		full.fit.uncertainty.fitParameterCovariance =
		    orbital_core::FitParameterCovariance::Identity(count, count) * 1e-8;
		QVERIFY(session.solution()->publish(full, "non-grav"));
		QCOMPARE(model->dimension(), count);
		QTest::qWait(120);
		auto* matrix = findItem(window->contentItem(), "covariance-matrix");
		QVERIFY(matrix);
		QVERIFY(matrix->property("contentWidth").toDouble() > matrix->width());
		model->selectCell(count - 1, count - 1);
		QCOMPARE(model->selection()["rowLabel"].toString(), QString("DT"));
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/covariance-nongrav.png"));
		model->setResult({});
		QTest::qWait(60);
		QVERIFY(model->selection().isEmpty());
		window->close();
	}

	void chartSnapshotBinsAndSelection()
	{
		auto result = fixture();
		result.observations.all.resize(10);
		result.fit.allResiduals.clear();
		result.fit.residuals.clear();
		result.fit.keptObsIndices = {0, 1, 2, 3, 4, 6, 7, 8, 9};
		const std::array values{0.0, 0.3, 0.5, 0.8, 1.2, 2.0, 0.1, 0.1, std::numeric_limits<double>::quiet_NaN(), 0.1};
		for (std::size_t i = 0; i < values.size(); ++i) {
			result.observations.all[i].JD = 2460000.5 + i;
			Residual r{};
			r.jd = 2461000.5 + i;
			r.dRA_arcsec = values[i];
			r.dDEC_arcsec = 0;
			r.total_arcsec = values[i];
			if (i == 6)
				result.observations.all[i].radar.emplace();
			if (i == 7)
				r.radar.emplace();
			if (i == 9)
				result.observations.all[i].JD = std::numeric_limits<double>::quiet_NaN();
			result.fit.allResiduals.add(r);
			if (i != 5)
				result.fit.residuals.add(r);
		}
		DiagnosticsModel model;
		model.setFromResult(result, 9);
		QScatterSeries ra, dec, excludedRa, excludedDec;
		auto stats = model.fillResidualCharts(&ra, &dec, &excludedRa, &excludedDec);
		QCOMPARE(stats["count"].toInt(), 6);
		QCOMPARE(stats["bins"].toList(), (QVariantList{1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0}));
		QCOMPARE(stats["excludedBins"].toList(), (QVariantList{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}));
		QCOMPARE(stats["histogramMax"].toDouble(), 2.0);
		QCOMPARE(stats["timeOriginMjd"].toDouble(), 60000.0);
		QCOMPARE(stats["timeOriginLabel"].toString(), QString("2023-02-25 UTC"));
		QCOMPARE(model.timeScale(), QString("UTC"));
		QCOMPARE(ra.at(0), QPointF(0, 0));
		QCOMPARE(excludedRa.at(0), QPointF(5, 2));
		QCOMPARE(dec.count(), 5);
		QCOMPARE(excludedDec.count(), 1);
		model.setFromResult(fixture(), 2);
		stats = model.fillResidualCharts(&ra, &dec, &excludedRa, &excludedDec);
		QCOMPARE(ra.count(), 2);
		QCOMPARE(excludedRa.count(), 1);
		QCOMPARE(stats["count"].toInt(), 3);
		auto partial = fixture();
		partial.fit.allResiduals.clear();
		model.setFromResult(partial, 2);
		stats = model.fillResidualCharts(&ra, &dec, &excludedRa, &excludedDec);
		QCOMPARE(stats["count"].toInt(), 2);
		QCOMPARE(ra.at(1).x(), 2.0);
		QVERIFY(excludedRa.points().isEmpty());
		model.clear();
		stats = model.fillResidualCharts(&ra, &dec, &excludedRa, &excludedDec);
		QCOMPARE(stats["count"].toInt(), 0);
		QCOMPARE(stats["bins"].toList(), (QVariantList{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
		QVERIFY(stats["jdMax"].toDouble() > stats["jdMin"].toDouble());
		QVERIFY(stats["draMax"].toDouble() > 0);
		QVERIFY(ra.points().isEmpty());
		QVERIFY(dec.points().isEmpty());
	}

	void largeResultPresentation()
	{
		if (!qEnvironmentVariableIsSet("ORBITAL_UI_BULK"))
			QSKIP("Large-result presentation benchmark is opt-in");
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window && QTest::qWaitForWindowExposed(window));
		auto result = fixture();
		result.observations.all.clear();
		result.fit.allResiduals.clear();
		result.fit.residuals.clear();
		result.fit.keptObsIndices.clear();
		for (int i = 0; i < 12000; ++i) {
			ObsNode observation;
			observation.JD = 2460000.5 + i * 0.04;
			observation.obsCode = i % 2 ? "A01" : "B02";
			result.observations.all.push_back(observation);
			Residual residual{};
			residual.dRA_arcsec = std::sin(i * 0.9) * 0.4;
			residual.dDEC_arcsec = std::cos(i * 1.1) * 0.3;
			residual.total_arcsec = std::hypot(residual.dRA_arcsec, residual.dDEC_arcsec);
			result.fit.allResiduals.add(residual);
			if (i % 13 != 0) {
				result.fit.keptObsIndices.push_back(i);
				result.fit.residuals.add(residual);
			}
		}
		QElapsedTimer timer;
		timer.start();
		QVERIFY(session.solution()->publish(result, "test"));
		const auto publishMs = timer.elapsed();
		QCoreApplication::processEvents();
		QVERIFY(!window->grabWindow().isNull());
		qInfo() << "UI_BULK observations=12000 publishMs=" << publishMs << "firstFrameMs=" << timer.elapsed();
		auto* ra = findItem(window->contentItem(), "ra-time-chart");
		QVERIFY(ra);
		auto* series = qobject_cast<QXYSeries*>(ra->property("includedSeries").value<QObject*>());
		QVERIFY(series && series->chart() && series->chart()->scene());
		const auto points = series->points();
		qInfo() << "UI_BULK chartSceneItems=" << series->chart()->scene()->items().size();
		QVERIFY(series->chart()->scene()->items().size() < 200);
		for (const QString theme : {"observatory", "light", "graphite"}) {
			timer.restart();
			session.setTheme(theme);
			const auto synchronousMs = timer.elapsed();
			QCoreApplication::processEvents();
			QVERIFY(!window->grabWindow().isNull());
			qInfo() << "UI_BULK theme=" << theme << "synchronousMs=" << synchronousMs << "frameMs=" << timer.elapsed();
			QCOMPARE(series->points(), points);
		}
		window->close();
	}

	void chartsLayoutAndInteraction()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		auto result = fixture();
		result.observations.all.clear();
		result.fit.allResiduals.clear();
		result.fit.residuals.clear();
		result.fit.keptObsIndices.clear();
		for (int i = 0; i < 240; ++i) {
			ObsNode o;
			o.JD = 2460000.5 + i * 0.4;
			o.obsCode = i % 2 ? "A01" : "B02";
			result.observations.all.push_back(o);
			Residual r{};
			r.dRA_arcsec = std::sin(i * 0.9) * (i % 13 == 0 ? 1.1 : 0.4);
			r.dDEC_arcsec = std::cos(i * 1.1) * (i % 13 == 0 ? 0.8 : 0.3);
			r.total_arcsec = std::hypot(r.dRA_arcsec, r.dDEC_arcsec);
			result.fit.allResiduals.add(r);
			if (i % 13 != 0) {
				result.fit.keptObsIndices.push_back(i);
				result.fit.residuals.add(r);
			}
		}
		QVERIFY(session.solution()->publish(result, "test"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QTest::qWait(100);
		auto* tab = findItem(window->contentItem(), "solution-charts-tab");
		QVERIFY(tab);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  tab->mapToScene({tab->width() / 2, tab->height() / 2}).toPoint());
		QCOMPARE(window->findChild<QObject*>("solution-page")->property("currentTab").toInt(), 2);
		auto* panel = findItem(window->contentItem(), "solution-charts");
		auto* ra = findItem(window->contentItem(), "ra-time-chart");
		auto* dec = findItem(window->contentItem(), "dec-time-chart");
		auto* hist = findItem(window->contentItem(), "residual-histogram");
		QVERIFY(panel && ra && dec && hist);
		auto* included = qobject_cast<QXYSeries*>(ra->property("includedSeries").value<QObject*>());
		auto* excluded = qobject_cast<QXYSeries*>(ra->property("excludedSeries").value<QObject*>());
		QVERIFY(included && excluded);
		QCOMPARE(included->count(), 221);
		QCOMPARE(excluded->count(), 19);
		QVERIFY(included->pointsVisible());
		QVERIFY(excluded->pointsVisible());
		QVERIFY(included->lightMarker().isNull());
		QVERIFY(excluded->lightMarker().isNull());
		window->setPosition(0, 0);
		window->requestActivate();
		QVERIFY(QTest::qWaitForWindowActive(window));
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		for (const QString theme : {"graphite", "observatory", "light"}) {
			session.setTheme(theme);
			for (const QSize size : {QSize(1024, 720), QSize(1512, 982), QSize(2554, 1312)}) {
				window->resize(size);
				panel->setProperty("contentY", 0);
				QTest::qWait(800);
				window->grabWindow();
				QTest::qWait(100);
				QVERIFY(panel->isVisible());
				QCOMPARE(included->chart()->size().toSize(), ra->size().toSize());
				QCOMPARE(panel->property("contentWidth").toDouble(), panel->width());
				const auto a = ra->mapRectToScene(ra->boundingRect());
				const auto d = dec->mapRectToScene(dec->boundingRect());
				const auto h = hist->mapRectToScene(hist->boundingRect());
				if (panel->property("stacked").toBool())
					QVERIFY(d.top() > a.bottom());
				else {
					QVERIFY(d.left() > a.right());
					QCOMPARE(a.top(), d.top());
				}
				QCOMPARE(ra->height(), dec->height());
				QVERIFY(ra->height() >= 230);
				QVERIFY(hist->height() >= 230);
				QVERIFY(h.top() > d.bottom());
				const auto contentHeight = panel->property("contentHeight").toDouble();
				QVERIFY(contentHeight + 1 >= panel->height());
				if (window->height() >= 982) {
					QVERIFY(std::abs(contentHeight - panel->height()) <= 1);
					QVERIFY(h.bottom() < panel->mapToScene({0, panel->height()}).y());
				} else if (window->height() <= 720) {
					QVERIFY(contentHeight > panel->height());
				}
				if (window->height() >= 1312) {
					QVERIFY(ra->height() > 320);
					QVERIFY(hist->height() > 300);
				}
				QVERIFY(d.right() <= window->width());
				QVERIFY(hist->width() > 0);
				const auto image = window->grabWindow();
				QVERIFY(!image.isNull());
				if (!captures.isEmpty())
					QVERIFY(image.save(captures + "/charts-" + theme + "-" + QString::number(size.width()) + "x" +
					                   QString::number(size.height()) + ".png"));
				const auto marker = ra->mapToScene(included->chart()->mapToPosition(excluded->at(1), excluded));
				const double scale = image.width() / static_cast<double>(window->width());
				const QPoint pixel(qRound(marker.x() * scale), qRound(marker.y() * scale));
				const QColor color = excluded->color();
				const QColor background = ra->property("backgroundColor").value<QColor>();
				bool markerVisible = false;
				int blendedEdgePixels = 0;
				const int radius = qCeil(4 * scale);
				for (int y = pixel.y() - radius; y <= pixel.y() + radius; ++y)
					for (int x = pixel.x() - radius; x <= pixel.x() + radius; ++x) {
						if (!image.rect().contains(x, y))
							continue;
						const auto actual = image.pixelColor(x, y);
						markerVisible |= std::abs(actual.red() - color.red()) < 20 &&
						                 std::abs(actual.green() - color.green()) < 20 &&
						                 std::abs(actual.blue() - color.blue()) < 20;
						const double coverage =
						    (actual.redF() - background.redF()) / (color.redF() - background.redF());
						if (coverage > 0.08 && coverage < 0.92 &&
						    std::abs(actual.greenF() - (background.greenF() +
						                                coverage * (color.greenF() - background.greenF()))) < 0.025 &&
						    std::abs(actual.blueF() -
						             (background.blueF() + coverage * (color.blueF() - background.blueF()))) < 0.025)
							++blendedEdgePixels;
					}
				QVERIFY2(markerVisible,
				         qPrintable(QString("Chart raster must follow point coordinates: %1, %2 px, marker (%3, %4)")
				                        .arg(theme)
				                        .arg(size.width())
				                        .arg(pixel.x())
				                        .arg(pixel.y())));
				QVERIFY2(blendedEdgePixels >= 4, "Residual markers must have antialiased edges");
			}
		}
		window->resize(1512, 982);
		window->grabWindow();
		QTest::qWait(100);
		window->requestActivate();
		QVERIFY(QTest::qWaitForWindowActive(window));
		QTest::mouseMove(window, QPoint(10, 10));
		QTRY_COMPARE(included->chart()->size().toSize(), ra->size().toSize());
		auto* observations = session.solution()->observations();
		const auto before = included->points();
		observations->setIncluded(1, false);
		observations->setFilterText("no match");
		QCOMPARE(included->points(), before);
		const auto point = excluded->at(1);
		const QPointF position = excluded->chart()->mapToPosition(point, excluded);
		ResidualMarker hitTest;
		QCOMPARE(hitTest.pointIndexAt(excluded, position), 1);
		QTest::mouseMove(window, ra->mapToScene(position).toPoint());
		auto* hint = findItem(ra, "residual-point-hint");
		QVERIFY(hint);
		QTRY_VERIFY(hint->isVisible());
		QCOMPARE(hint->property("observation").toMap()["sourceRow"].toInt(), 13);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier, ra->mapToScene(position).toPoint());
		QCOMPARE(observations->selectedSourceRow(), 13);
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/charts-hover.png"));
		QVERIFY(session.solution()->publish(fixture(), "next"));
		QTRY_COMPARE(included->count(), 2);
		QCOMPARE(excluded->count(), 1);
		session.solution()->diagnostics()->clear();
		QTRY_VERIFY(findItem(panel, "charts-empty")->isVisible());
		QCOMPARE(included->count(), 0);
		QCOMPARE(excluded->count(), 0);
		if (!captures.isEmpty())
			QVERIFY(window->grabWindow().save(captures + "/charts-empty.png"));
		window->close();
	}

	void residualViewsAndSourceSelection()
	{
		ObservationTableModel model;
		model.setResult(fixture());
		QScatterSeries included, excluded, selected;
		model.fillSeries(&included, true, false, false);
		model.fillSeries(&excluded, false, false, false);
		QCOMPARE(included.count(), 2);
		QCOMPARE(excluded.count(), 1);
		QCOMPARE(excluded.at(0), QPointF(0.1, -0.1));
		QCOMPARE(model.pointData(0.1, -0.1, false, false, false)["sourceRow"].toInt(), 1);
		QCOMPARE(model.selectedSourceRow(), -1);
		QVERIFY(model.plotBounds(false).contains(included.at(1)));
		model.selectPoint(0.1, -0.1, false, false);
		QCOMPARE(model.selectedSourceRow(), 1);
		QCOMPARE(model.selection()["station"].toString(), QString("B02"));
		model.setFilterText("A01");
		model.sortColumn(1, true);
		model.fillSeries(&excluded, false, true, false);
		QCOMPARE(excluded.at(0), QPointF(60001, 0.1));
		QCOMPARE(model.pointData(60001, 0.1, true, false, false)["station"].toString(), QString("B02"));
		QCOMPARE(model.pointData(60002, -0.2, true, true, true)["sourceRow"].toInt(), 2);
		QVERIFY(model.pointData(std::numeric_limits<double>::quiet_NaN(), 0, false, false, true).isEmpty());
		model.selectPoint(60002, -0.2, true, true);
		QCOMPARE(model.selectedSourceRow(), 2);
		model.fillSeries(&selected, true, false, false, true);
		QCOMPARE(selected.at(0), QPointF(0.2, -0.2));
		model.selectSourceRow(-1);
		QCOMPARE(model.selectedSourceRow(), 2);
		model.setResult({});
		QCOMPARE(model.stationResiduals()->count(), 0);
		QVERIFY(model.selection().isEmpty());
		QVERIFY(model.plotBounds(true).width() > 0);
		model.fillSeries(&selected, true, false, false, true);
		QCOMPARE(selected.count(), 0);
		QVERIFY(model.pointData(0, 0, false, false, true).isEmpty());
		auto overlap = fixture();
		overlap.observations.all[1].JD = overlap.observations.all[0].JD;
		overlap.fit.allResiduals = {};
		overlap.fit.allResiduals.add(Residual{});
		overlap.fit.allResiduals.add(Residual{});
		overlap.fit.allResiduals.add(Residual{});
		model.setResult(overlap);
		QCOMPARE(model.pointData(0, 0, false, false, false)["sourceRow"].toInt(), 1);
		QCOMPARE(model.pointData(0, 0, false, false, true)["sourceRow"].toInt(), 0);
	}

	void stationResidualStatistics()
	{
		StationResidualModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		const double nan = std::numeric_limits<double>::quiet_NaN();
		const std::vector<StationResidualSample> samples{
		    {"A01", true, QPointF(3, 4)},  {"A01", true, QPointF(-1, 2)},  {"A01", false, QPointF(100, 100)},
		    {"B02", false, QPointF(5, 5)}, {"C03", true, QPointF(nan, 1)}, {"C03", true, QPointF(0, 0)},
		    {"D04", true, std::nullopt}};
		model.setObservations(samples);
		QCOMPARE(model.count(), 4);
		QCOMPARE(model.included(), 5);
		QCOMPARE(model.excluded(), 2);
		QCOMPARE(model.missing(), 2);
		QCOMPARE(model.meanLimit(), 5.0);
		const auto value = [&](int row, int role) { return model.data(model.index(row, 0), role); };
		QCOMPARE(value(0, StationResidualModel::StationRole).toString(), QString("A01"));
		QCOMPARE(value(0, StationResidualModel::MeanRaRole).toDouble(), 1.0);
		QCOMPARE(value(0, StationResidualModel::MeanDecRole).toDouble(), 3.0);
		QVERIFY(std::abs(value(0, StationResidualModel::RmsRole).toDouble() - std::sqrt(7.5)) < 1e-14);
		QCOMPARE(value(1, StationResidualModel::ValidRole).toInt(), 1);
		QCOMPARE(value(1, StationResidualModel::RmsRole).toDouble(), 0.0);
		model.setSortMode(StationResidualModel::Rms);
		QCOMPARE(value(0, StationResidualModel::StationRole).toString(), QString("A01"));
		QVERIFY(std::isnan(value(3, StationResidualModel::RmsRole).toDouble()));
		model.setSortMode(StationResidualModel::MeanOffset);
		QCOMPARE(value(0, StationResidualModel::StationRole).toString(), QString("A01"));
		model.setSortMode(StationResidualModel::StationCode);
		QCOMPARE(value(1, StationResidualModel::StationRole).toString(), QString("B02"));
		model.setObservations({});
		QCOMPARE(model.count(), 0);
		QCOMPARE(model.missing(), 0);
		QCOMPARE(model.meanLimit(), 0.1);

		ObservationTableModel observations;
		auto result = fixture();
		result.fit.allResiduals = {};
		observations.setResult(result);
		auto* stations = observations.stationResiduals();
		QCOMPARE(stations->count(), 2);
		QCOMPARE(stations->included(), 2);
		QCOMPARE(stations->excluded(), 1);
		QCOMPARE(stations->missing(), 0);
		QCOMPARE(stations->data(stations->index(0, 0), StationResidualModel::MeanRaRole).toDouble(), 0.1);
		observations.setFilterText("2023");
		observations.setExcludedOnly(true);
		observations.sortColumn(1, true);
		QCOMPARE(stations->included(), 2);
		observations.showStation("A01");
		QCOMPARE(observations.count(), 2);
		QVERIFY(observations.filterText().isEmpty());
		QVERIFY(!observations.excludedOnly());
		QCOMPARE(observations.selectedSourceRow(), -1);
		observations.selectRow(0);
		QCOMPARE(observations.selection()["station"].toString(), QString("A01"));
		observations.setStationFilter("A0");
		QCOMPARE(observations.count(), 0);
		observations.showStation("B02");
		QCOMPARE(observations.count(), 1);
		observations.showStation("");
		QCOMPARE(observations.count(), 3);
		observations.setStationFilter("B02");
		result.observations.all[1].obsCode = "C03";
		observations.setResult(result);
		QVERIFY(observations.stationFilter().isEmpty());
		result.fit.residuals = {};
		observations.setResult(result);
		QCOMPARE(stations->missing(), 2);

		SolutionViewModel solution;
		QVERIFY(solution.publish(fixture(), "test"));
		auto invalid = fixture();
		invalid.fit.cancelled = true;
		QVERIFY(!solution.publish(invalid, "failed"));
		QCOMPARE(solution.observations()->stationResiduals()->included(), 2);
	}

	void normalizedRmsUsesAcceptedSolution()
	{
		SolutionViewModel solution;
		QVERIFY(std::isnan(solution.normalizedRms()));
		auto result = fixture();
		result.fit.normalizedRms = 1.25;
		QVERIFY(solution.publish(result, "test"));
		QCOMPARE(solution.rms(), result.fit.cost);
		QCOMPARE(solution.normalizedRms(), 1.25);
		result.fit.cancelled = true;
		result.fit.normalizedRms = 9.0;
		QVERIFY(!solution.publish(result, "cancelled"));
		QCOMPARE(solution.normalizedRms(), 1.25);
		QVERIFY(solution.publish(fixture(), "unavailable"));
		QVERIFY(std::isnan(solution.normalizedRms()));
	}

	void stationDiagnosticsAndCompactNavigation()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QVERIFY(session.solution()->publish(fixture(), "test"));
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		window->resize(1512, 982);
		QTRY_VERIFY(window->isVisible());
		const auto click = [&](QQuickItem* item) {
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
			                  item->mapToScene(QPointF(item->width() / 2, item->height() / 2)).toPoint());
		};
		auto* rail = findItem(window->contentItem(), "navigation-rail");
		auto* toggle = findItem(window->contentItem(), "toggle-navigation");
		auto* compactOpen = findItem(window->contentItem(), "compact-open-observations");
		auto* overview = findItem(window->contentItem(), "solution-overview");
		QVERIFY(rail && toggle && compactOpen && overview);
		auto* opticalMetric = findItem(window->contentItem(), "solution-optical-rms");
		auto* normalizedMetric = findItem(window->contentItem(), "solution-normalized-rms");
		QVERIFY(opticalMetric && normalizedMetric);
		QTRY_COMPARE(rail->width(), 202.0);
		click(toggle);
		QTRY_COMPARE(rail->width(), 59.0);
		QVERIFY(compactOpen->isVisible());
		auto* observationsButton = findItem(window->contentItem(), "navigation-page-1");
		QVERIFY(observationsButton);
		observationsButton->setProperty("count", "2777");
		click(observationsButton);
		QCOMPARE(window->property("currentPage").toInt(), 1);
		QVERIFY(observationsButton->property("checked").toBool());
		click(compactOpen);
		auto* importDialog = window->findChild<QObject*>("import-observations-dialog");
		QVERIFY(importDialog);
		QTRY_VERIFY(importDialog->property("visible").toBool());
		QVERIFY(QMetaObject::invokeMethod(importDialog, "close"));
		QTRY_VERIFY(!importDialog->property("visible").toBool());
		click(findItem(window->contentItem(), "navigation-page-0"));
		click(toggle);
		QTRY_COMPARE(rail->width(), 202.0);
		window->resize(1024, 720);
		QTRY_COMPARE(rail->width(), 59.0);
		QVERIFY(!toggle->isVisible());
		window->resize(1512, 982);
		QTRY_COMPARE(rail->width(), 202.0);
		const auto captures = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captures.isEmpty())
			QVERIFY(QDir().mkpath(captures));
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			for (const QSize size : {QSize(1024, 720), QSize(1512, 982), QSize(2554, 1312)}) {
				window->resize(size);
				QTest::qWait(80);
				QVERIFY(normalizedMetric->x() >= opticalMetric->x() + opticalMetric->width());
				QVERIFY(normalizedMetric->x() + normalizedMetric->width() <= normalizedMetric->parentItem()->width());
				QCOMPARE(normalizedMetric->y(), opticalMetric->y());
				for (auto* text : normalizedMetric->childItems())
					QVERIFY(!text->property("truncated").toBool());
				overview->setProperty(
				    "contentY", std::max(0.0, overview->property("contentHeight").toDouble() - overview->height()));
				QTest::qWait(60);
				auto* section = findItem(window->contentItem(), "station-diagnostics");
				QVERIFY(section);
				auto* navigationLabel = findItem(observationsButton, "navigation-label");
				auto* navigationCount = findItem(observationsButton, "navigation-count");
				QVERIFY(navigationLabel && navigationCount);
				if (!rail->property("compact").toBool()) {
					QVERIFY(!navigationLabel->property("truncated").toBool());
					QVERIFY(navigationLabel->x() + navigationLabel->width() < navigationCount->x());
				}
				const auto bounds = section->mapRectToScene(section->boundingRect());
				QVERIFY(bounds.left() >= rail->width());
				QVERIFY(bounds.right() <= window->width());
				QVERIFY(bounds.bottom() <= window->height());
				auto* row = findItem(window->contentItem(), "station-residual-row-A01");
				QVERIFY(row);
				QCOMPARE(row->property("includedCount").toInt(), 2);
				auto* rms = findItem(row, "station-rms-value");
				auto* rmsHeader = findItem(section, "station-rms-header");
				QVERIFY(rms && rmsHeader);
				QVERIFY(!rms->property("truncated").toBool());
				const auto rmsBounds = rms->mapRectToScene(rms->boundingRect());
				const auto headerBounds = rmsHeader->mapRectToScene(rmsHeader->boundingRect());
				QVERIFY(std::abs(rmsBounds.right() - headerBounds.right()) < 0.1);
				const auto rowBounds = row->mapRectToScene(row->boundingRect());
				QVERIFY(std::abs(rmsBounds.right() - (rowBounds.right() - 14)) < 0.1);
				int offsetCount = 0;
				const auto inspect = [&](const auto& self, QQuickItem* item) -> void {
					if (item->objectName() == "residual-offset-value") {
						++offsetCount;
						QVERIFY(!item->property("truncated").toBool());
						QVERIFY(item->mapRectToScene(item->boundingRect()).right() < rmsBounds.left());
					}
					for (auto* child : item->childItems())
						self(self, child);
				};
				inspect(inspect, row);
				QCOMPARE(offsetCount, 2);
				if (!captures.isEmpty())
					QVERIFY(window->grabWindow().save(captures + "/stations-" + theme + "-" +
					                                  QString::number(size.width()) + ".png"));
			}
		}
		window->resize(1512, 982);
		QTest::qWait(80);
		overview->setProperty("contentY",
		                      std::max(0.0, overview->property("contentHeight").toDouble() - overview->height()));
		QTest::qWait(60);
		auto* stationLink = findItem(window->contentItem(), "station-link-A01");
		QVERIFY(stationLink && stationLink->isVisible());
		session.solution()->observations()->setFilterText("no match");
		session.solution()->observations()->setExcludedOnly(true);
		click(stationLink);
		QTRY_COMPARE(window->property("currentPage").toInt(), 1);
		QCOMPARE(session.solution()->observations()->stationFilter(), QString("A01"));
		QCOMPARE(session.solution()->observations()->count(), 2);
		auto* stationFilter = findItem(window->contentItem(), "observations-station-filter");
		QVERIFY(stationFilter);
		QCOMPARE(stationFilter->property("currentIndex").toInt(), 1);
		window->setProperty("currentPage", 0);
		QTest::qWait(50);
		auto* sort = findItem(window->contentItem(), "station-diagnostics-sort");
		QVERIFY(sort);
		sort->setProperty("currentIndex", int(StationResidualModel::Rms));
		QVERIFY(QMetaObject::invokeMethod(sort, "activated", Q_ARG(int, int(StationResidualModel::Rms))));
		QCOMPARE(session.solution()->observations()->stationResiduals()->sortMode(), StationResidualModel::Rms);
		auto* all = findItem(window->contentItem(), "all-station-observations");
		QVERIFY(all);
		click(all);
		QCOMPARE(window->property("currentPage").toInt(), 1);
		QCOMPARE(session.solution()->observations()->count(), 3);
		window->close();
	}

	void controlCursors_data()
	{
		QTest::addColumn<QByteArray>("declaration");
		QTest::addColumn<int>("cursor");
		QTest::newRow("button") << QByteArray("CustomButton { text: 'Run' }") << int(Qt::PointingHandCursor);
		QTest::newRow("checkbox") << QByteArray("CustomCheckBox { text: 'Include' }") << int(Qt::PointingHandCursor);
		QTest::newRow("tab") << QByteArray("CustomTab { text: 'Overview' }") << int(Qt::PointingHandCursor);
		QTest::newRow("link") << QByteArray("CustomTextLink { text: 'Open table' }") << int(Qt::PointingHandCursor);
		QTest::newRow("navigation") << QByteArray("NavigationItem { text: 'Solution' }") << int(Qt::PointingHandCursor);
		QTest::newRow("segments") << QByteArray("CustomSegmentedControl { options: ['Scatter', 'By time'] }")
		                          << int(Qt::PointingHandCursor);
		QTest::newRow("compact-navigation")
		    << QByteArray("NavigationItem { text: 'Solution'; compact: true }") << int(Qt::PointingHandCursor);
		QTest::newRow("input") << QByteArray("CustomTextField {}") << int(Qt::IBeamCursor);
		QTest::newRow("scrollbar")
		    << QByteArray("CustomScrollBar { orientation: Qt.Horizontal; size: 0.5; policy: CustomScrollBar.AlwaysOn }")
		    << int(Qt::PointingHandCursor);
		QTest::newRow("static-scrollbar")
		    << QByteArray("CustomScrollBar { orientation: Qt.Horizontal; size: 1; policy: CustomScrollBar.AlwaysOn }")
		    << int(Qt::ArrowCursor);
	}

	void controlCursors()
	{
		QFETCH(QByteArray, declaration);
		QFETCH(int, cursor);
		QQmlEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		QQmlComponent component(&engine);
		component.setData("import QtQuick\nimport TrajectaUI\n" + declaration, QUrl());
		QScopedPointer<QObject> object(component.create());
		QVERIFY2(object, qPrintable(component.errorString()));
		auto* item = qobject_cast<QQuickItem*>(object.data());
		QVERIFY(item);
		QQuickWindow window;
		window.resize(240, 100);
		item->setParentItem(window.contentItem());
		item->setPosition({20, 20});
		item->setSize({180, 40});
		window.show();
		QVERIFY(QTest::qWaitForWindowExposed(&window));
		QTest::mouseMove(&window, {5, 5});
		QTest::mouseMove(&window, {60, 40});
		QTRY_COMPARE(int(window.cursor().shape()), cursor);
		QTest::mouseMove(&window, {5, 5});
		QTRY_COMPARE(window.cursor().shape(), Qt::ArrowCursor);
		item->setEnabled(false);
		QTest::mouseMove(&window, {60, 40});
		QTRY_COMPARE(window.cursor().shape(), Qt::ArrowCursor);
		item->setParentItem(nullptr);
	}

	void weightingDefaultsAndSavedValues()
	{
		QTemporaryDir dir;
		QVERIFY(dir.isValid());
		QSettings store(dir.filePath("settings.ini"), QSettings::IniFormat);
		WeightingSettings weighting;
		SigmaClipSettings clipping;
		weighting.load(store);
		clipping.load(store);
		QVERIFY(weighting.useWeights());
		QVERIFY(weighting.observerSigmas());
		QVERIFY(weighting.debiasEnabled());
		QVERIFY(!weighting.timingEnabled());
		QCOMPARE(weighting.weightProfile(), 1);
		QVERIFY(clipping.enabled());
		QCOMPARE(clipping.threshold(), 3.0);
		store.setValue("useWeights", false);
		store.setValue("observerSigmas", false);
		store.setValue("debiasEnabled", false);
		store.setValue("timingEnabled", true);
		store.setValue("weightProfile", 2);
		store.setValue("sigmaClipEnabled", false);
		store.setValue("sigmaClipThreshold", 4.0);
		weighting.load(store);
		clipping.load(store);
		QVERIFY(!weighting.useWeights());
		QVERIFY(!weighting.observerSigmas());
		QVERIFY(!weighting.debiasEnabled());
		QVERIFY(weighting.timingEnabled());
		QCOMPARE(weighting.weightProfile(), 2);
		QVERIFY(!clipping.enabled());
		QCOMPARE(clipping.threshold(), 4.0);
	}

	void calculationSettingsTransaction()
	{
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		QSettings().clear();
		QSettings().setValue("forces/sunJ2", false);
		FitSettingsAggregator settings;
		settings.load();
		QVERIFY(!QSettings().contains("forces/sunJ2"));
		ApplicationSession session(settings);
		settings.weighting()->setProperty("useWeights", false);
		settings.perturbations()->setProperty("enabled", false);
		settings.save();
		auto* draft = session.calculationSettings();
		draft->reset();
		QVERIFY(draft->valid());
		draft->weighting()->setProperty("useWeights", true);
		draft->sigmaClip()->setProperty("threshold", 4.5);
		QVERIFY(!settings.weighting()->useWeights());
		QCOMPARE(settings.sigmaClip()->threshold(), 3.0);
		QVERIFY(!QSettings().value("fit/useWeights").toBool());
		session.setTheme("slate");
		draft->reset();
		QVERIFY(!draft->weighting()->useWeights());
		QCOMPARE(session.theme(), QString("slate"));
		QCOMPARE(QSettings().value("appearance/theme").toString(), QString("slate"));
		draft->weighting()->setProperty("useWeights", true);
		draft->sigmaClip()->setProperty("threshold", 4.5);
		draft->perturbations()->setProperty("integratorInitialStepDays", 2.5);
		draft->perturbations()->setProperty("relativityMode", 1);
		draft->perturbations()->setProperty("earthJ3", false);
		draft->perturbations()->setProperty("earthJ4", false);
		draft->perturbations()->setProperty("rotatingEarthPole", false);
		draft->photometry()->setProperty("asteroidG", 0.42);
		draft->photometry()->setProperty("cometK1", 8.0);
		draft->photometry()->setProperty("cometK2", 5.0);
		draft->photometry()->setProperty("nuclearPhaseCoefficient", 0.03);
		draft->nonGrav()->setProperty("model", 4);
		QCOMPARE(draft->nonGravDefaultMask(), 2);
		QCOMPARE(draft->nonGravSupportedMask(), 2);
		draft->nonGrav()->setProperty("solveMask", 2);
		QVERIFY(!session.applySettings(true));
		QVERIFY(!settings.weighting()->useWeights());
		QSignalSpy changes(&settings, &FitSettingsAggregator::anyChanged);
		QVERIFY(session.applySettings());
		QCOMPARE(changes.count(), 1);
		QVERIFY(settings.weighting()->useWeights());
		QCOMPARE(settings.sigmaClip()->threshold(), 4.5);
		QCOMPARE(QSettings().value("integrator/initialStepDays").toDouble(), 2.5);
		QCOMPARE(QSettings().value("forces/relativityMode").toInt(), 1);
		QCOMPARE(QSettings().value("non-gravitational/model").toInt(), 4);
		QCOMPARE(QSettings().value("non-gravitational/solveMask").toInt(), 2);
		QCOMPARE(QSettings().value("photometry/asteroidG").toDouble(), 0.42);
		QCOMPARE(QSettings().value("photometry/cometK1").toDouble(), 8.0);
		QCOMPARE(QSettings().value("photometry/cometK2").toDouble(), 5.0);
		QCOMPARE(QSettings().value("photometry/nuclearPhaseCoefficient").toDouble(), 0.03);
		FitSettingsAggregator restored;
		restored.load();
		const auto restoredForces = restored.perturbations()->forceModel();
		QVERIFY(restoredForces.relativity == orbital_core::RelativityMode::SunOuter);
		QVERIFY(restoredForces.earthJ2);
		QVERIFY(!restoredForces.earthJ3);
		QVERIFY(!restoredForces.earthJ4);
		QVERIFY(!restoredForces.rotatingEarthPole);
		QCOMPARE(restored.photometry()->options().asteroidG, 0.42);
		QCOMPARE(restored.photometry()->options().cometK1, 8.0);
		QCOMPARE(restored.photometry()->options().cometK2, 5.0);
		QCOMPARE(restored.photometry()->options().cometNuclearPhaseCoefficient, 0.03);
		for (double value :
		     {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
			draft->sigmaClip()->setProperty("threshold", value);
			QVERIFY(!draft->valid());
			QVERIFY(!session.applySettings());
			QCOMPARE(settings.sigmaClip()->threshold(), 4.5);
			draft->reset();
		}
		draft->nonGrav()->setProperty("solveMask", 1);
		QVERIFY(!draft->valid());
		draft->reset();
		draft->nonGrav()->setProperty("model", 1);
		draft->nonGrav()->setProperty("solveMask", 0);
		QCOMPARE(draft->nonGravDefaultMask(), 15);
		draft->nonGrav()->setProperty("solveMask", 2);
		QVERIFY(!draft->valid());
		draft->nonGrav()->setProperty("solveMask", 3);
		QVERIFY(draft->valid());
		draft->perturbations()->setProperty("integratorInitialStepDays", 1e-8);
		QVERIFY(!draft->valid());
	}

	void observatoryUpdateChecksBeforeWriting()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const auto previous = orbital_core::observatoryDataDirectory();
		const auto restore = qScopeGuard([&] { orbital_core::setObservatoryDataDirectory(previous); });
		orbital_core::setObservatoryDataDirectory(directory.path().toStdWString());
		QFile file(directory.filePath("OBSCODE.dat"));
		QVERIFY(file.open(QIODevice::WriteOnly));
		QCOMPARE(file.write("original"), 8);
		file.close();
		ObservatoryUpdate updater([](const QString& id, std::string_view current, std::stop_token) {
			if (id == "radar")
				throw std::runtime_error("offline");
			return orbital_core::obs::ObservatoryUpdate{"updated", 1, 0, current != "updated"};
		});
		QSignalSpy installed(&updater, &ObservatoryUpdate::installed);
		updater.check();
		QTRY_VERIFY(!updater.states().value("observatories").toMap().value("checking").toBool());
		QTRY_VERIFY(!updater.states().value("radar").toMap().value("checking").toBool());
		QCOMPARE(updater.availableCount(), 1);
		QVERIFY(!updater.states().value("radar").toMap().value("message").toString().isEmpty());
		QVERIFY(file.open(QIODevice::ReadOnly));
		QCOMPARE(file.readAll(), QByteArray("original"));
		file.close();
		updater.install("observatories");
		QCOMPARE(installed.size(), 1);
		QCOMPARE(updater.availableCount(), 0);
		QVERIFY(file.open(QIODevice::ReadOnly));
		QCOMPARE(file.readAll(), QByteArray("updated"));
		file.close();
		updater.check();
		QTRY_VERIFY(!updater.states().value("observatories").toMap().value("checking").toBool());
		QTRY_VERIFY(!updater.states().value("radar").toMap().value("checking").toBool());
		QCOMPARE(updater.availableCount(), 0);
		QCOMPARE(installed.size(), 1);
	}

	void earthPckUpdateChecksAndRejectsInvalidFile()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QTcpServer server;
		QVERIFY(server.listen(QHostAddress::LocalHost));
		QByteArray payload = QByteArray("DAF/PCK") + QByteArray(2048, 'x');
		const QByteArray listing =
		    "earth_620120_260806.bpc earth_000101_270208_261112.bpc "
		    "earth_2026_260806_2126_predict.bpc";
		connect(&server, &QTcpServer::newConnection, &server, [&] {
			auto* socket = server.nextPendingConnection();
			connect(socket, &QTcpSocket::readyRead, socket, [socket, &listing, &payload] {
				const auto request = socket->readAll();
				const auto& body = request.startsWith("GET / ") ? listing : payload;
				socket->write("HTTP/1.1 200 OK\r\nContent-Length: " + QByteArray::number(body.size()) +
				              "\r\nConnection: close\r\n\r\n" + body);
				socket->disconnectFromHost();
			});
		});
		EarthPckUpdate updater(QUrl(QString("http://127.0.0.1:%1/").arg(server.serverPort())), directory.path(),
		                       directory.filePath("updates.ini"));
		QSignalSpy installed(&updater, &EarthPckUpdate::installed);
		updater.check();
		QTRY_VERIFY(!updater.states().value("earth_pck_latest").toMap().value("checking").toBool());
		QCOMPARE(updater.availableCount(), 3);
		updater.install("earth_pck_latest");
		QTRY_VERIFY(!updater.states().value("earth_pck_latest").toMap().value("downloading").toBool());
		QCOMPARE(installed.size(), 0);
		QVERIFY(!QFile::exists(directory.filePath("earth_latest.bpc")));
		QVERIFY(!updater.states().value("earth_pck_latest").toMap().value("message").toString().isEmpty());

		QFile valid(QString::fromStdWString(orbital_core::earthPckFilePath("earth_latest.bpc").wstring()));
		if (!valid.open(QIODevice::ReadOnly))
			QSKIP("Earth PCK fixture unavailable");
		payload = valid.readAll();
		updater.install("earth_pck_latest");
		QTRY_VERIFY_WITH_TIMEOUT(!updater.states().value("earth_pck_latest").toMap().value("downloading").toBool(),
		                         30000);
		QCOMPARE(installed.size(), 1);
		QVERIFY(!QFile::exists(directory.filePath("earth_latest.bpc")));
		QVERIFY(QFile::exists(directory.filePath("earth_latest.bpc.pending")));
		QFile changedFile(directory.filePath("earth_latest.bpc"));
		QVERIFY(changedFile.open(QIODevice::WriteOnly));
		QCOMPARE(changedFile.write("changed"), 7);
		changedFile.close();
		EarthPckUpdate::applyPending(directory.path(), directory.filePath("updates.ini"));
		QVERIFY(changedFile.open(QIODevice::ReadOnly));
		QCOMPARE(changedFile.readAll(), QByteArray("changed"));
		changedFile.close();
		QVERIFY(changedFile.remove());
		EarthPckUpdate::applyPending(directory.path(), directory.filePath("updates.ini"));
		QFile output(directory.filePath("earth_latest.bpc"));
		QVERIFY(output.open(QIODevice::ReadOnly));
		QCOMPARE(output.readAll(), payload);
		QVERIFY(!QFile::exists(directory.filePath("earth_latest.bpc.pending")));
		updater.check();
		QTRY_VERIFY(!updater.states().value("earth_pck_latest").toMap().value("checking").toBool());
		QCOMPARE(updater.availableCount(), 2);
	}

	void earthPckUpdateKeepsOpenFileUntilRestart()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QFile fixture(QString::fromStdWString(orbital_core::earthPckFilePath("earth_latest.bpc").wstring()));
		if (!fixture.open(QIODevice::ReadOnly))
			QSKIP("Earth PCK fixture unavailable");
		const QByteArray payload = fixture.readAll();
		const auto local = directory.filePath("earth_latest.bpc");
		QFile old(local);
		QVERIFY(old.open(QIODevice::WriteOnly));
		QCOMPARE(old.write(payload), payload.size());
		QCOMPARE(old.write("x", 1), 1);
		old.close();
		QVERIFY(old.open(QIODevice::ReadOnly));
		QTcpServer server;
		QVERIFY(server.listen(QHostAddress::LocalHost));
		const QByteArray listing =
		    "earth_620120_260806.bpc earth_000101_270209_261113.bpc "
		    "earth_2026_260806_2126_predict.bpc";
		connect(&server, &QTcpServer::newConnection, &server, [&] {
			auto* socket = server.nextPendingConnection();
			connect(socket, &QTcpSocket::readyRead, socket, [socket, &listing, &payload] {
				const auto request = socket->readAll();
				const auto& body = request.startsWith("GET / ") ? listing : payload;
				socket->write("HTTP/1.1 200 OK\r\nContent-Length: " + QByteArray::number(body.size()) +
				              "\r\nConnection: close\r\n\r\n" + body);
				socket->disconnectFromHost();
			});
		});
		const auto metadata = directory.filePath("updates.ini");
		EarthPckUpdate updater(QUrl(QString("http://127.0.0.1:%1/").arg(server.serverPort())), directory.path(),
		                       metadata);
		QSignalSpy installed(&updater, &EarthPckUpdate::installed);
		updater.check();
		QTRY_VERIFY(!updater.states().value("earth_pck_latest").toMap().value("checking").toBool());
		QVERIFY(updater.states().value("earth_pck_latest").toMap().value("available").toBool());
		updater.install("earth_pck_latest");
		QTRY_VERIFY_WITH_TIMEOUT(!updater.states().value("earth_pck_latest").toMap().value("downloading").toBool(),
		                         30000);
		QCOMPARE(installed.size(), 1);
		QVERIFY(QFile::exists(local + ".pending"));
		QCOMPARE(old.size(), payload.size() + 1);
		old.close();
		EarthPckUpdate::applyPending(directory.path(), metadata);
		QVERIFY(old.open(QIODevice::ReadOnly));
		QCOMPARE(old.readAll(), payload);
		QVERIFY(!QFile::exists(local + ".pending"));
	}

	void earthPckUpdatedFileTakesPriority()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const auto previousBase = getExecutableDir();
		const auto previousPck = orbital_core::earthPckDataDirectory();
		const auto restore = qScopeGuard([&] {
			getExecutableDir() = previousBase;
			orbital_core::setEarthPckDataDirectory(previousPck);
		});
		getExecutableDir() = directory.filePath("app").toStdWString();
		orbital_core::setEarthPckDataDirectory(directory.filePath("user").toStdWString());
		QVERIFY(QDir().mkpath(directory.filePath("app/earth_orientation")));
		QVERIFY(QDir().mkpath(directory.filePath("user")));
		QFile bundled(directory.filePath("app/earth_orientation/earth_latest.bpc"));
		QVERIFY(bundled.open(QIODevice::WriteOnly));
		bundled.close();
		QCOMPARE(QDir::fromNativeSeparators(
		             QString::fromStdWString(orbital_core::earthPckFilePath("earth_latest.bpc").wstring())),
		         directory.filePath("user/earth_latest.bpc"));
		QCOMPARE(QDir::fromNativeSeparators(QString::fromStdWString(orbital_core::earthC04DataDirectory().wstring())),
		         directory.filePath("user"));
		QFile updated(directory.filePath("user/earth_latest.bpc"));
		QVERIFY(updated.open(QIODevice::WriteOnly));
		updated.close();
		QCOMPARE(QFileInfo(QString::fromStdWString(orbital_core::earthPckFilePath("earth_latest.bpc").wstring()))
		             .canonicalFilePath(),
		         QFileInfo(updated).canonicalFilePath());
	}

	void observatoryUpdateDoesNotOverwriteConcurrentChange()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		const auto previous = orbital_core::observatoryDataDirectory();
		const auto restore = qScopeGuard([&] { orbital_core::setObservatoryDataDirectory(previous); });
		orbital_core::setObservatoryDataDirectory(directory.path().toStdWString());
		QFile file(directory.filePath("OBSCODE.dat"));
		QVERIFY(file.open(QIODevice::WriteOnly));
		QCOMPARE(file.write("original"), 8);
		file.close();
		ObservatoryUpdate updater([](const QString& id, std::string_view, std::stop_token) {
			return orbital_core::obs::ObservatoryUpdate{"updated", 1, 0, id == "observatories"};
		});
		updater.check();
		QTRY_VERIFY(!updater.states().value("observatories").toMap().value("checking").toBool());
		QTRY_VERIFY(!updater.states().value("radar").toMap().value("checking").toBool());
		QVERIFY(file.open(QIODevice::WriteOnly | QIODevice::Truncate));
		QCOMPARE(file.write("edited"), 6);
		file.close();
		updater.install("observatories");
		QVERIFY(file.open(QIODevice::ReadOnly));
		QCOMPARE(file.readAll(), QByteArray("edited"));
		file.close();
		QTRY_VERIFY(!updater.states().value("observatories").toMap().value("checking").toBool());
		QCOMPARE(updater.availableCount(), 1);
	}

	void observatoryUpdateControls()
	{
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QTRY_VERIFY(window->isVisible());
		auto* button = findItem(window->contentItem(), "data-sources-button");
		QVERIFY(button);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  button->mapToScene({button->width() / 2, button->height() / 2}).toPoint());
		QTRY_VERIFY(session.dataSources()->checked());
		QTRY_VERIFY(findItem(window->contentItem(), "update-source-observatories"));
		QVERIFY(findItem(window->contentItem(), "update-source-radar"));
		QVERIFY(findItem(window->contentItem(), "update-source-earth_pck_historic"));
		QVERIFY(findItem(window->contentItem(), "update-source-earth_pck_latest"));
		QVERIFY(findItem(window->contentItem(), "update-source-earth_pck_predict"));
	}

	void calculationSettingsDialog()
	{
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		QSettings().clear();
		QSettings().setValue("integrator/nBodyPropagation", false);
		FitSettingsAggregator settings;
		settings.weighting()->setProperty("debiasEnabled", false);
		ApplicationSession session(settings);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QTRY_VERIFY(window->isVisible());
		auto* dialog = window->findChild<QObject*>("calculation-settings-dialog");
		auto* backdrop = window->findChild<QObject*>("dialog-backdrop");
		QVERIFY(dialog && backdrop);
		const auto click = [&](const QString& name) {
			auto* item = findItem(window->contentItem(), name);
			if (!item || !item->isVisible() || !item->isEnabled())
				return false;
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
			                  item->mapToScene({item->width() / 2, item->height() / 2}).toPoint());
			return true;
		};
		QVERIFY(click("settings-button"));
		QTRY_VERIFY(dialog->property("visible").toBool());
		QVERIFY(backdrop->property("active").toBool());
		QTest::qWait(100);
		QVERIFY(click("theme-slate"));
		QCOMPARE(session.theme(), QString("slate"));
		QVERIFY(click("settings-section-1"));
		QTest::qWait(100);
		auto* weightProfile = findItem(window->contentItem(), "settings-weight-profile");
		QVERIFY(weightProfile);
		QCOMPARE(weightProfile->property("count").toInt(), 2);
		QCOMPARE(weightProfile->property("currentIndex").toInt(), 0);
		session.calculationSettings()->weighting()->setProperty("useWeights", true);
		QVERIFY(click("settings-weight-profile"));
		QTest::keyClick(window, Qt::Key_End);
		QTest::keyClick(window, Qt::Key_Return);
		QCOMPARE(session.calculationSettings()->weighting()->weightProfile(), 2);
		QCOMPARE(weightProfile->property("currentText").toString(), QString("Uniform 0.5 arcsec"));
		QVERIFY(click("settings-weight-profile"));
		QTest::keyClick(window, Qt::Key_Home);
		QTest::keyClick(window, Qt::Key_Return);
		QCOMPARE(session.calculationSettings()->weighting()->weightProfile(), 1);
		QVERIFY(click("settings-debias"));
		QVERIFY(session.calculationSettings()->weighting()->debiasEnabled());
		QVERIFY(!settings.weighting()->debiasEnabled());
		QVERIFY(click("settings-cancel"));
		QTRY_VERIFY(!dialog->property("visible").toBool());
		QVERIFY(!settings.weighting()->debiasEnabled());
		QCOMPARE(session.theme(), QString("slate"));
		QVERIFY(click("settings-button"));
		QTest::qWait(100);
		QVERIFY(!session.calculationSettings()->weighting()->debiasEnabled());
		QVERIFY(click("settings-debias"));
		QVERIFY(click("settings-apply"));
		QTRY_VERIFY(!dialog->property("visible").toBool());
		QVERIFY(settings.weighting()->debiasEnabled());
		QVERIFY(click("settings-button"));
		QTest::qWait(100);
		if (!session.calculationSettings()->sigmaClip()->enabled())
			QVERIFY(click("settings-clipping"));
		auto* threshold = findItem(window->contentItem(), "settings-clip-threshold");
		QVERIFY(threshold);
		QVERIFY(click("settings-clip-threshold"));
		QTest::keyClick(window, Qt::Key_A, Qt::ControlModifier);
		QTest::keyClick(window, Qt::Key_Minus);
		QVERIFY(!session.calculationSettings()->valid());
		QVERIFY(!findItem(window->contentItem(), "settings-apply")->isEnabled());
		QTest::keyClick(window, Qt::Key_A, Qt::ControlModifier);
		QTest::keyClick(window, Qt::Key_4);
		QVERIFY(session.calculationSettings()->valid());
		QCOMPARE(session.calculationSettings()->sigmaClip()->threshold(), 4.0);
		QTest::keyClick(window, Qt::Key_Escape);
		QTRY_VERIFY(!dialog->property("visible").toBool());
		QCOMPARE(settings.sigmaClip()->threshold(), 3.0);
		QVERIFY(click("settings-button"));
		QTest::qWait(100);
		QVERIFY(click("settings-section-4"));
		QTest::qWait(100);
		QVERIFY(click("settings-nongrav"));
		QVERIFY(session.calculationSettings()->nonGrav()->enabled());
		QVERIFY(click("settings-nongrav-model"));
		QTest::keyClick(window, Qt::Key_End);
		QTest::keyClick(window, Qt::Key_Return);
		QCOMPARE(session.calculationSettings()->nonGrav()->model(), 4);
		QCOMPARE(session.calculationSettings()->nonGravDefaultMask(), 2);
		QVERIFY(click("settings-nongrav-model"));
		QTest::keyClick(window, Qt::Key_Up);
		QTest::keyClick(window, Qt::Key_Return);
		QCOMPARE(session.calculationSettings()->nonGrav()->model(), 3);
		QCOMPARE(session.calculationSettings()->nonGravDefaultMask(), 7);
		QVERIFY(click("settings-nongrav-coefficients"));
		QTest::keyClick(window, Qt::Key_End);
		QTest::keyClick(window, Qt::Key_Return);
		QCOMPARE(session.calculationSettings()->nonGrav()->solveMask(), 15);
		QVERIFY(click("settings-cancel"));
		const auto captureDir = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captureDir.isEmpty())
			QVERIFY(QDir().mkpath(captureDir));
		for (const auto size : {QSize(1512, 982), QSize(1024, 720)}) {
			window->resize(size);
			QTest::qWait(150);
			session.setTheme(size.width() > 1024 ? "graphite" : "light");
			QVERIFY(click("settings-button"));
			for (int section = 0; section < 6; ++section) {
				dialog->setProperty("currentSection", section);
				QTest::qWait(150);
				QVERIFY(dialog->property("height").toDouble() <= window->height() - 40);
				if (section == 4) {
					auto* model = findItem(window->contentItem(), "settings-nongrav-model");
					QVERIFY(model && model->width() >= 120);
					QVERIFY(!findItem(window->contentItem(), "a1-prior"));
					QVERIFY(!findItem(window->contentItem(), "a1-sigma"));
				}
				auto* apply = findItem(window->contentItem(), "settings-apply");
				QVERIFY(apply);
				const auto point = apply->mapToScene({apply->width(), apply->height()});
				QVERIFY(point.x() < window->width() && point.y() < window->height());
				if (!captureDir.isEmpty())
					QVERIFY(window->grabWindow().save(captureDir + "/settings-" + QString::number(section) + "-" +
					                                  QString::number(size.width()) + ".png"));
			}
			QVERIFY(click("settings-cancel"));
		}
		window->close();
	}

	void nativeShellAndLiveFit()
	{
		const auto oldExecutableDir = getExecutableDir();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source de440;
		SpkReader sb441;
		setDefaultEphemerisSource(&de440);
		setDefaultAsteroidSource(&sb441);
		const auto restoreSources = qScopeGuard([&] {
			setDefaultEphemerisSource(nullptr);
			setDefaultAsteroidSource(nullptr);
			getExecutableDir() = oldExecutableDir;
		});
		FitSettingsAggregator settings;
		const auto inputUrl =
		    QUrl::fromLocalFile(QDir::cleanPath(QStringLiteral(ORBITAL_UI_SOURCE_DIR) +
		                                        "/../OrbitalCore/tests/data/2020VV_ADES_2022.ESA11_table_replay.xml"));
		QSettings().setValue("recentFiles", QStringList{inputUrl.toLocalFile(), m_settings.filePath("missing.xml")});
		ApplicationSession session(settings);
		QCOMPARE(session.recentFiles()->count(), 1);
		QCOMPARE(session.recentFiles()->filePath(0), inputUrl.toLocalFile());
		QSignalSpy runs(&session, &ApplicationSession::runStarted);
		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QTRY_VERIFY(window->isVisible());
		QVERIFY(!session.canRun());
		QVERIFY(!session.openInput(QUrl("https://example.invalid/observations.xml")));
		QVERIFY(!session.error().isEmpty());
		QVERIFY(!session.canRun());
		session.dismissError();
		const auto captureDir = qEnvironmentVariable("ORBITAL_UI_CAPTURES");
		if (!captureDir.isEmpty())
			QVERIFY(QDir().mkpath(captureDir));
		auto* dropArea = window->findChild<QQuickItem*>("observation-drop-area");
		QVERIFY(dropArea);
		QTRY_VERIFY(dropArea->isEnabled());
		const auto dropInput = [&](const QList<QUrl>& urls, QPoint position) {
			QMimeData mime;
			mime.setUrls(urls);
			QDragEnterEvent enter(position, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
			QCoreApplication::sendEvent(window, &enter);
			if (!enter.isAccepted()) {
				QDragLeaveEvent leave;
				QCoreApplication::sendEvent(window, &leave);
				return false;
			}
			QDragMoveEvent move(position, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
			QCoreApplication::sendEvent(window, &move);
			QDropEvent drop(position, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
			QCoreApplication::sendEvent(window, &drop);
			return drop.isAccepted() && drop.dropAction() == Qt::CopyAction;
		};
		const QList<QList<QUrl>> rejectedDrops{
		    {},
		    {QUrl("https://example.invalid/observations.xml")},
		    {inputUrl, inputUrl},
		    {QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR))},
		    {QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/missing.xml")}};
		for (const auto& urls : rejectedDrops) {
			QVERIFY(!session.canDropInput(urls));
			QVERIFY(!dropInput(urls, {700, 450}));
		}
		QVERIFY(!session.canRun());
		QVERIFY(session.canDropInput({inputUrl}));
		for (const int pageIndex : {0, 1}) {
			window->setProperty("currentPage", pageIndex);
			for (const QPoint position : {QPoint(30, 120), QPoint(420, 30), QPoint(700, 450), QPoint(1000, 720)}) {
				QMimeData mime;
				mime.setUrls({inputUrl});
				QDragEnterEvent enter(position, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
				QCoreApplication::sendEvent(window, &enter);
				QVERIFY(enter.isAccepted());
				QDragLeaveEvent leave;
				QCoreApplication::sendEvent(window, &leave);
				QVERIFY(!session.canRun());
				QVERIFY(!session.busy());
				QVERIFY(!session.solution()->available());
				QVERIFY(!dropArea->property("containsDrag").toBool());
			}
		}
		window->setProperty("currentPage", 0);
		const auto selectedInput = session.inputFile();
		for (const auto& urls : rejectedDrops) {
			QVERIFY(!dropInput(urls, {700, 450}));
			QCOMPARE(session.inputFile(), selectedInput);
		}
		QMimeData dragMime;
		dragMime.setUrls({inputUrl});
		QDragEnterEvent dragEnter({700, 450}, Qt::CopyAction, &dragMime, Qt::LeftButton, Qt::NoModifier);
		QCoreApplication::sendEvent(window, &dragEnter);
		QVERIFY(dragEnter.isAccepted());
		QVERIFY(dropArea->property("containsDrag").toBool());
		QTest::qWait(50);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/observations-drop.png"));
		QDragLeaveEvent dragLeave;
		QCoreApplication::sendEvent(window, &dragLeave);
		QVERIFY(!dropArea->property("containsDrag").toBool());
		auto* importDialog = window->findChild<QObject*>("import-observations-dialog");
		auto* calculation = window->findChild<QObject*>("calculation-dialog");
		auto* backdrop = window->findChild<QObject*>("dialog-backdrop");
		QVERIFY(importDialog && calculation && backdrop);
		auto* openButton = findItem(window->contentItem(), "open-observations-button");
		QVERIFY(openButton);
		window->requestActivate();
		QTRY_COMPARE(QGuiApplication::focusWindow(), window);
		QTest::mouseMove(window, {window->width() - 20, 60});
		QTest::qWait(50);
		QTest::mouseMove(window, openButton->mapToScene({20, 20}).toPoint());
		QTRY_VERIFY(openButton->property("hovered").toBool());
		QTest::qWait(800);
		QVERIFY(!backdrop->property("active").toBool());
		QVERIFY(dropArea->isEnabled());
		QVERIFY(QMetaObject::invokeMethod(importDialog, "open"));
		QTest::qWait(100);
		QVERIFY(!dropArea->isEnabled());
		QVERIFY(backdrop->property("active").toBool());
		auto* recentRow = findItem(window->contentItem(), "recent-observation-0");
		QVERIFY(recentRow);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/import-observations.png"));
		auto* picker = window->findChild<QObject*>("observations-file-picker");
		auto* chooseButton = findItem(window->contentItem(), "choose-observations-button");
		QVERIFY(picker && chooseButton);
		QVERIFY(!QCoreApplication::testAttribute(Qt::AA_DontUseNativeDialogs));
		QVERIFY(!(picker->property("options").toInt() & QFileDialog::DontUseNativeDialog));
		QVERIFY(!picker->property("visible").toBool());
		QVERIFY(!session.busy());
		QVERIFY(QMetaObject::invokeMethod(importDialog, "close"));
		QTRY_VERIFY(dropArea->isEnabled());
		if (qEnvironmentVariableIsSet("ORBITAL_UI_LIVE_FIT")) {
			settings.perturbations()->setProperty("enabled", true);
			settings.perturbations()->setUseSb441(true);
			settings.weighting()->setProperty("useWeights", true);
			settings.weighting()->setProperty("debiasEnabled", true);
			settings.weighting()->setProperty("observerSigmas", false);
			settings.weighting()->setProperty("timingEnabled", false);
			settings.sigmaClip()->setProperty("enabled", true);
			QVERIFY(dropInput({inputUrl}, {700, 450}));
			QVERIFY(session.busy());
			QCOMPARE(runs.count(), 1);
			QVERIFY(calculation->property("visible").toBool());
			QVERIFY(!calculation->property("closeAllowed").toBool());
			QCOMPARE(QFileInfo(session.inputFile()).canonicalFilePath(),
			         QFileInfo(inputUrl.toLocalFile()).canonicalFilePath());
			RecentFiles recentReader;
			QCOMPARE(recentReader.count(), 1);
			QCOMPARE(recentReader.filePath(0), session.inputFile());
			QVERIFY(!session.canDropInput({inputUrl}));
			QVERIFY(!dropInput({inputUrl}, {700, 450}));
			session.run();
			QCOMPARE(runs.count(), 1);
			QTest::qWait(40);
			if (!captureDir.isEmpty())
				QVERIFY(window->grabWindow().save(captureDir + "/calculation.png"));
			QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
			QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
			QVERIFY(session.solution()->available());
			QVERIFY(session.solution()->included() > 3);
			QVERIFY(std::isfinite(session.solution()->absoluteMagnitude()));
			QVERIFY(std::isfinite(session.solution()->earthMoidKm()));
			const auto& orbit = *session.solution()->snapshot()->fit.orbit;
			const auto moids = orbital_core::hazard::computeMoidsToAllPlanets(orbit, orbit.getEpoch());
			QCOMPARE(session.solution()->earthMoidKm(),
			         moids[orbital_core::hazard::Planet::Earth] * (Constants::AU_IN_METERS / 1000.0));
			QVERIFY(QDir(session.solution()->directory()).exists());
			QVERIFY(QFile::exists(session.solution()->directory() + "/orbit_cartesian_oef.txt"));
			QVERIFY(QFile::exists(session.solution()->directory() + "/weighted_residuals_rwo.txt"));
			QVERIFY(QMetaObject::invokeMethod(calculation, "close"));
			QTRY_VERIFY(dropArea->isEnabled());
		} else {
			QVERIFY(session.solution()->publish(fixture(), "fixture"));
		}
		session.solution()->observations()->selectSourceRow(0);
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			auto* runButton = window->findChild<QQuickItem*>("run-button");
			QVERIFY(runButton);
			const bool enabled = runButton->isEnabled();
			const QColor expected = enabled ? QColor(theme == "light" ? "#ffffff" : "#06120f")
			                                : QColor(theme == "light" ? "#64736f" : "#8d9998");
			QCOMPARE(runButton->property("foreground").value<QColor>(), expected);
			for (const QSize size :
			     {QSize(1024, 720), QSize(1280, 800), QSize(1512, 982), QSize(1920, 1080), QSize(2554, 1312)}) {
				window->resize(size);
				QTest::qWait(120);
				const auto image = window->grabWindow();
				QVERIFY(!image.isNull());
				const auto scale = image.width() / static_cast<double>(window->width());
				QVERIFY(image.pixelColor(qRound(5 * scale), qRound(70 * scale)) !=
				        image.pixelColor(qRound(230 * scale), qRound(70 * scale)));
				if (!captureDir.isEmpty())
					QVERIFY(
					    image.save(captureDir + "/solution-" + theme + "-" + QString::number(size.width()) + ".png"));
			}
		}
		window->resize(1512, 982);
		session.setTheme("graphite");
		QTest::qWait(150);
		if (qEnvironmentVariableIsSet("ORBITAL_UI_3D")) {
			window->setProperty("currentPage", 2);
			QTest::qWait(200);
			auto* scene = findItem(window->contentItem(), "orbit-scene");
			QVERIFY(scene);
			QVERIFY(session.solution()->scene()->error().isEmpty());
			if (!captureDir.isEmpty())
				QVERIFY(window->grabWindow().save(captureDir + "/orbit-live-fit.png"));
			window->setProperty("currentPage", 0);
			QTest::qWait(80);
		}
		auto* elementsMode = window->findChild<QQuickItem*>("orbit-elements-mode");
		auto* equButton = findItem(window->contentItem(), "orbit-elements-mode-option-1");
		auto* kepButton = findItem(window->contentItem(), "orbit-elements-mode-option-0");
		QVERIFY(elementsMode && equButton && kepButton);
		QVERIFY(equButton->isEnabled());
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier, equButton->mapToScene({20, 12}).toPoint());
		QTRY_COMPARE(elementsMode->property("currentIndex").toInt(), 1);
		auto* secondElement = findItem(window->contentItem(), "orbit-element-1");
		QVERIFY(secondElement);
		QCOMPARE(secondElement->property("modelData").toMap()["label"].toString(), QString("h"));
		QTest::qWait(100);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/solution-equinoctial.png"));
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier, kepButton->mapToScene({20, 12}).toPoint());
		QTRY_COMPARE(elementsMode->property("currentIndex").toInt(), 0);
		secondElement = findItem(window->contentItem(), "orbit-element-1");
		QVERIFY(secondElement);
		QCOMPARE(secondElement->property("modelData").toMap()["label"].toString(), QString("e"));
		auto* plot = window->findChild<QQuickItem*>("residual-plot");
		QVERIFY(plot);
		QXYSeries* series = nullptr;
		for (int i = 0; i < plot->property("count").toInt(); ++i) {
			QAbstractSeries* candidate = nullptr;
			QVERIFY(
			    QMetaObject::invokeMethod(plot, "series", Q_RETURN_ARG(QAbstractSeries*, candidate), Q_ARG(int, i)));
			if (candidate && candidate->objectName() == "included-residuals")
				series = qobject_cast<QXYSeries*>(candidate);
		}
		QVERIFY(series && series->count() > 0);
		for (int i = 0; i < plot->property("count").toInt(); ++i) {
			QAbstractSeries* candidate = nullptr;
			QVERIFY(
			    QMetaObject::invokeMethod(plot, "series", Q_RETURN_ARG(QAbstractSeries*, candidate), Q_ARG(int, i)));
			if (!candidate || !candidate->objectName().startsWith("residual-"))
				continue;
			auto* circle = qobject_cast<QLineSeries*>(candidate);
			QVERIFY(circle);
			QCOMPARE(circle->count(), 129);
			const double multiplier = circle->objectName() == "residual-rms-circle" ? 1 : 3;
			const double radius = session.solution()->rms() * multiplier;
			QCOMPARE(circle->at(0), QPointF(radius, 0));
			QPointF horizontalPoint, verticalPoint, origin;
			QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, horizontalPoint),
			                                  Q_ARG(QPointF, QPointF(radius, 0)), Q_ARG(QAbstractSeries*, candidate)));
			QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, verticalPoint),
			                                  Q_ARG(QPointF, QPointF(0, radius)), Q_ARG(QAbstractSeries*, candidate)));
			QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, origin),
			                                  Q_ARG(QPointF, QPointF(0, 0)), Q_ARG(QAbstractSeries*, candidate)));
			QVERIFY(std::abs(QLineF(origin, horizontalPoint).length() - QLineF(origin, verticalPoint).length()) < 1.0);
		}
		const auto points = series->points();
		const auto point = *std::max_element(points.cbegin(), points.cend(),
		                                     [](const auto& a, const auto& b) { return a.x() < b.x(); });
		QPointF pointPosition;
		QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, pointPosition),
		                                  Q_ARG(QPointF, point), Q_ARG(QAbstractSeries*, series)));
		QTest::mouseMove(window, plot->mapToScene({1, 1}).toPoint());
		QTest::mouseMove(window, plot->mapToScene(pointPosition).toPoint());
		QTRY_VERIFY_WITH_TIMEOUT(series->property("pointHovered").toBool(), 1500);
		QTRY_COMPARE(window->cursor().shape(), Qt::PointingHandCursor);
		auto* pointHint = window->findChild<QQuickItem*>("residual-point-hint");
		QVERIFY(pointHint);
		QTRY_VERIFY(pointHint->isVisible());
		QCOMPARE(pointHint->property("observation").toMap(),
		         session.solution()->observations()->pointData(point.x(), point.y(), false, false, true));
		QVERIFY(!backdrop->property("active").toBool());
		QVERIFY(pointHint->x() >= 0 && pointHint->x() + pointHint->width() <= plot->width());
		QVERIFY(pointHint->y() >= 0 && pointHint->y() + pointHint->height() <= plot->height());
		QTest::qWait(100);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/residuals-hover.png"));
		QTest::mouseMove(window, plot->mapToScene({-8, -8}).toPoint());
		QTRY_COMPARE(window->cursor().shape(), Qt::ArrowCursor);
		QTRY_VERIFY(!series->property("pointHovered").toBool());
		QVERIFY(!pointHint->isVisible());
		QTest::mouseMove(window, plot->mapToScene({1, 1}).toPoint());
		QTRY_COMPARE(window->cursor().shape(), Qt::ArrowCursor);
		auto* mode = window->findChild<QObject*>("residual-mode");
		QVERIFY(mode);
		QVERIFY(QMetaObject::invokeMethod(mode, "activated", Q_ARG(int, 1)));
		QQmlListReference axes(plot, "axes");
		QCOMPARE(axes.count(), 2);
		QCOMPARE(plot->property("timeMode").toBool(), true);
		QTest::qWait(150);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/residuals-time.png"));
		const auto timePoint = series->at(series->count() - 1);
		QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, pointPosition),
		                                  Q_ARG(QPointF, timePoint), Q_ARG(QAbstractSeries*, series)));
		QTest::mouseMove(window, plot->mapToScene(pointPosition).toPoint());
		QTRY_VERIFY(pointHint->isVisible());
		QCOMPARE(pointHint->property("observation").toMap(),
		         session.solution()->observations()->pointData(timePoint.x(), timePoint.y(), true, false, true));
		QVERIFY(!backdrop->property("active").toBool());
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/residuals-time-hover.png"));
		QXYSeries* decSeries = nullptr;
		for (int i = 0; i < plot->property("count").toInt(); ++i) {
			QAbstractSeries* candidate = nullptr;
			QVERIFY(
			    QMetaObject::invokeMethod(plot, "series", Q_RETURN_ARG(QAbstractSeries*, candidate), Q_ARG(int, i)));
			if (candidate && candidate->objectName() == "included-dec-residuals")
				decSeries = qobject_cast<QXYSeries*>(candidate);
		}
		QVERIFY(decSeries && decSeries->count() > 0);
		const auto decPoint = decSeries->at(decSeries->count() - 1);
		QVERIFY(QMetaObject::invokeMethod(plot, "mapToPosition", Q_RETURN_ARG(QPointF, pointPosition),
		                                  Q_ARG(QPointF, decPoint), Q_ARG(QAbstractSeries*, decSeries)));
		QTest::mouseMove(window, plot->mapToScene(pointPosition).toPoint());
		QTRY_VERIFY(pointHint->isVisible());
		QCOMPARE(pointHint->property("observation").toMap(),
		         session.solution()->observations()->pointData(decPoint.x(), decPoint.y(), true, true, true));
		QVERIFY(QMetaObject::invokeMethod(mode, "activated", Q_ARG(int, 0)));
		QVERIFY(!pointHint->isVisible());
		auto* overview = window->findChild<QQuickItem*>("solution-overview");
		QVERIFY(overview);
		auto* stationDiagnostics = findItem(window->contentItem(), "station-diagnostics");
		QVERIFY(stationDiagnostics);
		overview->setProperty("contentY",
		                      std::max(0.0, overview->property("contentHeight").toDouble() - overview->height()));
		QTest::qWait(150);
		session.solution()->observations()->selectSourceRow(1);
		QCOMPARE(session.solution()->observations()->selectedSourceRow(), 1);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/overview-scrolled.png"));
		auto* observationLink = findItem(window->contentItem(), "observation-details-link");
		auto* observationDetails = window->findChild<QObject*>("observation-details-dialog");
		QVERIFY(observationLink && observationDetails);
		overview->setProperty("contentY", std::max(0.0, observationLink->mapToItem(overview, QPointF()).y() +
		                                                    overview->property("contentY").toDouble() - 40));
		QTest::qWait(50);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  observationLink->mapToScene({20, observationLink->height() / 2}).toPoint());
		QTRY_VERIFY(observationDetails->property("visible").toBool());
		QCOMPARE(observationDetails->property("title").toString(), QString("Observation 2"));
		QCOMPARE(window->property("currentPage").toInt(), 0);
		QVERIFY(backdrop->property("active").toBool());
		QTest::qWait(100);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/observation-details.png"));
		window->resize(1024, 720);
		QTest::qWait(100);
		QVERIFY(observationDetails->property("height").toDouble() <= window->height() - 40);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/observation-details-compact.png"));
		QVERIFY(QMetaObject::invokeMethod(observationDetails, "close"));
		QTRY_VERIFY(!backdrop->property("active").toBool());
		auto* moidDialog = window->findChild<QObject*>("moid-details-dialog");
		QVERIFY(moidDialog);
		QVERIFY(QMetaObject::invokeMethod(moidDialog, "open"));
		QTest::qWait(100);
		QVERIFY(backdrop->property("active").toBool());
		QVERIFY(moidDialog->property("height").toDouble() < window->height() - 40);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/moid-details.png"));
		QVERIFY(QMetaObject::invokeMethod(moidDialog, "close"));
		QTRY_VERIFY(!backdrop->property("active").toBool());
		auto* details = window->findChild<QObject*>("orbit-details-dialog");
		QVERIFY(details);
		for (const QString theme : {"graphite", "light"}) {
			session.setTheme(theme);
			for (const QSize size : {QSize(1024, 720), QSize(1512, 982)}) {
				window->resize(size);
				QVERIFY(QMetaObject::invokeMethod(details, "open"));
				QTest::qWait(150);
				QVERIFY(details->property("width").toDouble() <= window->width() - 56);
				QVERIFY(details->property("height").toDouble() <= window->height() - 56);
				auto* content = details->property("contentItem").value<QQuickItem*>();
				QVERIFY(content);
				auto* derived = findItem(content, "derived-orbit-parameters");
				auto* state = findItem(content, "cartesian-state-table");
				QVERIFY(derived && state);
				QVERIFY(derived->x() + derived->width() <= state->x());
				QVERIFY(state->x() + state->width() <= content->width());
				auto* sizeSection = findItem(content, "estimated-size");
				QVERIFY(sizeSection);
				const bool hasMagnitude = std::isfinite(session.solution()->absoluteMagnitude());
				QCOMPARE(sizeSection->isVisible(), hasMagnitude);
				QVERIFY(sizeSection->x() + sizeSection->width() <= content->width());
				const auto sizeRows = session.solution()->sizeEstimates();
				QCOMPARE(sizeRows.size(), hasMagnitude ? 5 : 0);
				for (const auto& row : sizeRows) {
					auto* item = findItem(sizeSection, "size-row-" + row.toMap()["label"].toString());
					QVERIFY(item);
					auto* value = findItem(item, "size-value");
					QVERIFY(value);
					QCOMPARE(value->property("text").toString(), row.toMap()["value"].toString());
					QVERIFY(value->width() >= value->property("contentWidth").toDouble());
				}
				QVERIFY(findItem(sizeSection, "typical-size"));
				QVERIFY(!dropArea->isEnabled());
				QVERIFY(!dropInput({inputUrl}, {700, 450}));
				if (!captureDir.isEmpty())
					QVERIFY(window->grabWindow().save(captureDir + "/orbit-details-" + theme + "-" +
					                                  QString::number(size.width()) + ".png"));
				auto* header = details->property("header").value<QQuickItem*>();
				QVERIFY(header);
				auto* closeButton = findItem(header, "dialog-close-button");
				QVERIFY(closeButton);
				QTest::mouseClick(
				    window, Qt::LeftButton, Qt::NoModifier,
				    closeButton->mapToScene(QPointF(closeButton->width() / 2, closeButton->height() / 2)).toPoint());
				QTRY_VERIFY(!details->property("visible").toBool());
				QTRY_VERIFY(dropArea->isEnabled());
			}
		}
		window->resize(1512, 982);
		session.setTheme("graphite");
		auto* page = window->findChild<QObject*>("solution-page");
		QVERIFY(page);
		page->setProperty("currentTab", 1);
		QTest::qWait(150);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/covariance.png"));
		auto* appearance = window->findChild<QObject*>("calculation-settings-dialog");
		QVERIFY(appearance);
		QVERIFY(QMetaObject::invokeMethod(appearance, "open"));
		QTest::qWait(150);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/appearance.png"));
		QVERIFY(QMetaObject::invokeMethod(appearance, "close"));
		if (qEnvironmentVariableIsSet("ORBITAL_UI_LIVE_FIT")) {
			QVERIFY(QMetaObject::invokeMethod(appearance, "open"));
			QTest::qWait(100);
			auto* rerunButton = findItem(window->contentItem(), "settings-recalculate");
			QVERIFY(rerunButton && rerunButton->isEnabled());
			session.calculationSettings()->perturbations()->setProperty("integratorInitialStepDays", 5.0);
			const auto beforeRerun = runs.count();
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier, rerunButton->mapToScene({20, 16}).toPoint());
			QCOMPARE(runs.count(), beforeRerun + 1);
			QVERIFY(!session.applySettings());
			QVERIFY(!session.applySettings(true));
			QVERIFY(calculation->property("visible").toBool());
			QVERIFY(!appearance->property("visible").toBool());
			QCOMPARE(settings.perturbations()->integratorInitialStepDays(), 5.0);
			QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 120000);
			QVERIFY2(session.error().isEmpty(), qPrintable(session.error()));
			QVERIFY(session.solution()->available());
			QVERIFY(QMetaObject::invokeMethod(calculation, "close"));
		}
		window->setProperty("currentPage", 1);
		session.solution()->observations()->selectRow(0);
		QTest::qWait(150);
		auto* table = window->findChild<QQuickItem*>("observations-table");
		QVERIFY(table);
		QVERIFY(table->width() >= 700);
		if (!captureDir.isEmpty())
			QVERIFY(window->grabWindow().save(captureDir + "/observations.png"));
		auto* filter = window->findChild<QQuickItem*>("observations-filter");
		QVERIFY(filter);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  filter->mapToScene({20, filter->height() / 2}).toPoint());
		QTest::qWait(100);
		for (const char letter : QByteArray("no-matching-station"))
			QTest::keyClick(window, letter);
		QCOMPARE(session.solution()->observations()->count(), 0);
		QTest::keyClick(window, Qt::Key_A, Qt::ControlModifier);
		QTest::keyClick(window, Qt::Key_Backspace);
		QCOMPARE(session.solution()->observations()->count(), session.solution()->total());
		const auto previousTitle = session.solution()->objectTitle();
		const auto previousDirectory = session.solution()->directory();
		QTemporaryDir invalid;
		QFile invalidFile(invalid.filePath("broken.xml"));
		QVERIFY(invalidFile.open(QIODevice::WriteOnly));
		invalidFile.write("not observations");
		invalidFile.close();
		session.recentFiles()->addFile(QUrl::fromLocalFile(invalidFile.fileName()));
		QVERIFY(QMetaObject::invokeMethod(importDialog, "open"));
		QTest::qWait(50);
		recentRow = findItem(window->contentItem(), "recent-observation-0");
		QVERIFY(recentRow);
		const auto beforeRecent = runs.count();
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  recentRow->mapToScene({20, recentRow->height() / 2}).toPoint());
		QCOMPARE(runs.count(), beforeRecent + 1);
		QVERIFY(calculation->property("visible").toBool());
		QTRY_VERIFY(!importDialog->property("visible").toBool());
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 10000);
		QVERIFY(!session.error().isEmpty());
		QVERIFY(calculation->property("visible").toBool());
		QVERIFY(calculation->property("closeAllowed").toBool());
		QCOMPARE(session.solution()->objectTitle(), previousTitle);
		QCOMPARE(session.solution()->directory(), previousDirectory);
		QVERIFY(QMetaObject::invokeMethod(calculation, "close"));
		const auto previousRunCount = runs.count();
		QVERIFY(!session.openInput(QUrl::fromLocalFile(invalid.filePath("missing.xml"))));
		QCOMPARE(runs.count(), previousRunCount);
		QVERIFY(!calculation->property("visible").toBool());
		QVERIFY(!session.openRecent(-1));
		QVERIFY(!session.openRecent(100));
		QVERIFY(session.openInput(inputUrl));
		auto* footer = calculation->property("footer").value<QQuickItem*>();
		QVERIFY(footer);
		auto* cancelButton = findItem(footer, "dialog-accept-button");
		QVERIFY(cancelButton);
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier, cancelButton->mapToScene({20, 16}).toPoint());
		QVERIFY(!session.busy() || session.cancelling());
		QVERIFY(calculation->property("visible").toBool());
		QTRY_VERIFY_WITH_TIMEOUT(!session.busy(), 10000);
		QCOMPARE(session.progress(), QString("Cancelled"));
		QVERIFY(session.error().isEmpty());
		QCOMPARE(session.solution()->directory(), previousDirectory);
		QVERIFY(QMetaObject::invokeMethod(calculation, "close"));
		window->close();
	}

	void ephemerisColumnSelection()
	{
		EphemerisTableModel model;
		QAbstractItemModelTester tester(&model, QAbstractItemModelTester::FailureReportingMode::QtTest);
		QCOMPARE(model.columns(), 7);
		QVERIFY(model.includesColumn(EphemerisTableModel::Distance));
		QVERIFY(model.includesColumn(EphemerisTableModel::SolarDistance));
		QVERIFY(model.includesColumn(EphemerisTableModel::Magnitude));
		QVERIFY(model.includesColumn(EphemerisTableModel::Elongation));
		QVERIFY(!model.includesColumn(EphemerisTableModel::LongitudeRate));
		QVERIFY(!model.includesColumn(EphemerisTableModel::Phase));
		QVERIFY(!model.uncertaintySelected());
		QStringList optionKeys;
		for (const auto& option : model.columnOptions())
			optionKeys << option.toMap().value("key").toString();
		QVERIFY(optionKeys.contains("ellipseMajor"));
		QVERIFY(!optionKeys.contains("ellipseMinor"));
		QVERIFY(!optionKeys.contains("ellipseAngle"));
		QVERIFY(!optionKeys.contains("rangeUncertainty"));
		QVERIFY(!optionKeys.contains("longitudeUncertainty"));
		model.setColumnEnabled(EphemerisTableModel::EllipseMajor, true);
		QVERIFY(model.includesColumn(EphemerisTableModel::EllipseMajor));
		QVERIFY(model.includesColumn(EphemerisTableModel::EllipseMinor));
		QVERIFY(model.includesColumn(EphemerisTableModel::EllipseAngle));
		model.setColumnEnabled(EphemerisTableModel::EllipseMajor, false);
		QVERIFY(!model.uncertaintySelected());
		orbital_core::EphemerisRow row;
		row.mjd = 60000;
		row.longitudeRate = 3.0 / (206264.80624709636 / 24.0);
		row.latitudeRate = 4.0 / (206264.80624709636 / 24.0);
		model.setRows({row}, orbital_core::EphemerisPlane::Equatorial, orbital_core::time::Scale::UTC, false);
		QVERIFY(model.header(-1).isEmpty());
		QVERIFY(model.header(model.columns()).isEmpty());
		QVERIFY(model.cell(0, -1).isEmpty());
		QVERIFY(model.cell(0, model.columns()).isEmpty());
		model.setColumnEnabled(EphemerisTableModel::Distance, false);
		model.setColumnEnabled(EphemerisTableModel::SkyMotion, true);
		model.setColumnEnabled(EphemerisTableModel::SkyPositionAngle, true);
		QCOMPARE(model.columns(), 8);
		QVERIFY(!model.includesColumn(EphemerisTableModel::Distance));
		QCOMPARE(model.cell(0, 6), QString("5.000"));
		QCOMPARE(model.cell(0, 7), QString("36.870"));
		EphemerisTableModel reloaded;
		QVERIFY(reloaded.includesColumn(EphemerisTableModel::Distance));
		QVERIFY(!reloaded.includesColumn(EphemerisTableModel::SkyMotion));
		model.setLayout(EphemerisTableModel::Layout::Cartesian);
		QCOMPARE(model.columns(), 7);
		QVERIFY(!model.includesColumn(EphemerisTableModel::SkyMotion));
		model.setColumnEnabled(EphemerisTableModel::SkyMotion, true);
		QCOMPARE(model.columns(), 8);
		QCOMPARE(model.header(7), QString::fromUtf8("Sky motion, \u2033/h"));
		QCOMPARE(model.cell(0, 7), QString("5.000"));
		model.setLayout(EphemerisTableModel::Layout::Spherical);
		QCOMPARE(model.columns(), 8);
		row.longitudeRate = row.latitudeRate = 0.0;
		model.setRows({row}, orbital_core::EphemerisPlane::Equatorial, orbital_core::time::Scale::UTC, false);
		QCOMPARE(model.cell(0, 7), QString::fromUtf8("\u2014"));
	}

	void ephemeridesPage()
	{
		const auto oldExecutableDir = getExecutableDir();
		getExecutableDir() = QCoreApplication::applicationDirPath().toStdWString();
		DE440Source de440;
		QVERIFY2(de440.isAvailable(), "Test planetary ephemeris was not loaded");
		setDefaultEphemerisSource(&de440);
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(nullptr);
			getExecutableDir() = oldExecutableDir;
		});
		FitSettingsAggregator settings;
		ApplicationSession session(settings);
		const auto result = fixture();
		QVERIFY(session.solution()->publish(result, "test"));
		const QTemporaryDir exports;
		QVERIFY(exports.isValid());
		auto* ephemerides = session.ephemerides();
		QVERIFY(!ephemerides->available());
		ephemerides->setSolution(result, "Fixture", 18.0, nullptr, exports.path());
		QVERIFY(ephemerides->available());
		QAbstractItemModelTester tester(ephemerides->table(), QAbstractItemModelTester::FailureReportingMode::QtTest);
		QCOMPARE(ephemerides->startDate(), QString("2023-02-24"));
		QCOMPARE(ephemerides->endDate(), QString("2023-03-10"));
		QCOMPARE(ephemerides->plannedEpochs(), 15);
		QVERIFY(ephemerides->formError().isEmpty());

		QQmlApplicationEngine engine;
		engine.setImportPathList({QLibraryInfo::path(QLibraryInfo::QmlImportsPath)});
		engine.setInitialProperties({{"session", QVariant::fromValue(&session)}});
		engine.load(QUrl::fromLocalFile(QStringLiteral(ORBITAL_UI_SOURCE_DIR) + "/qml/shell/WorkspaceWindow.qml"));
		QVERIFY(!engine.rootObjects().isEmpty());
		auto* window = qobject_cast<QQuickWindow*>(engine.rootObjects().front());
		QVERIFY(window);
		QVERIFY(QTest::qWaitForWindowExposed(window));
		const auto click = [&](const QString& name) {
			auto* item = findItem(window->contentItem(), name);
			QVERIFY2(item, qPrintable(name));
			QVERIFY2(item->isVisible() && item->isEnabled(), qPrintable(name));
			QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
			                  item->mapToScene({item->width() / 2, item->height() / 2}).toPoint());
		};
		click("navigation-page-3");
		auto* page = findItem(window->contentItem(), "ephemerides-page");
		QVERIFY(page);
		QTRY_VERIFY(page->isVisible());
		QCOMPARE(findItem(page, "ephemeris-planned-epochs")->property("text").toString(), QString("15 epochs"));
		QVERIFY(findItem(page, "ephemeris-station")->isVisible());
		click("ephemeris-settings-button");
		auto* ephemerisDialog = window->findChild<QObject*>("ephemeris-settings-dialog");
		QVERIFY(ephemerisDialog);
		QTRY_VERIFY(ephemerisDialog->property("visible").toBool());
		auto* dialogContent = ephemerisDialog->property("contentItem").value<QQuickItem*>();
		QVERIFY(dialogContent);
		auto* motionOption = findItem(dialogContent, "ephemeris-column-skyMotion");
		QVERIFY(motionOption);
		QVERIFY(!motionOption->property("checked").toBool());
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  motionOption->mapToScene({motionOption->width() / 2, motionOption->height() / 2}).toPoint());
		QVERIFY(ephemerides->table()->includesColumn(EphemerisTableModel::SkyMotion));
		motionOption = findItem(dialogContent, "ephemeris-column-skyMotion");
		QVERIFY(motionOption);
		QTRY_VERIFY(motionOption->y() > 0);
		QVERIFY(motionOption->property("checked").toBool());
		QTest::mouseClick(window, Qt::LeftButton, Qt::NoModifier,
		                  motionOption->mapToScene({motionOption->width() / 2, motionOption->height() / 2}).toPoint());
		QVERIFY(!ephemerides->table()->includesColumn(EphemerisTableModel::SkyMotion));
		QVERIFY(QMetaObject::invokeMethod(ephemerisDialog, "close"));

		click("ephemeris-generate");
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QCOMPARE(ephemerides->epochCount(), 15);
		const auto firstArchive = ephemerides->archivePath();
		QVERIFY(firstArchive.startsWith(exports.filePath("Ephemerides/run-")));
		QVERIFY(QFileInfo::exists(QDir(firstArchive).filePath("ephemeris.txt")));
		QVERIFY(QFileInfo::exists(QDir(firstArchive).filePath("result.json")));
		QFile autoJson(QDir(firstArchive).filePath("result.json"));
		QVERIFY(autoJson.open(QIODevice::ReadOnly));
		const auto autoDocument = QJsonDocument::fromJson(autoJson.readAll()).object();
		QCOMPARE(autoDocument.value("rows").toArray().size(), 15);
		QVERIFY(autoDocument.value("rows").toArray().first().toObject().value("range_3sigma_km").isDouble());
		QFile autoText(QDir(firstArchive).filePath("ephemeris.txt"));
		QVERIFY(autoText.open(QIODevice::ReadOnly));
		QVERIFY(autoText.readAll().contains("$$SOE"));
		QVERIFY(findItem(page, "ephemeris-open-archive")->isEnabled());
		QVERIFY(!findItem(page, "ephemeris-export"));
		ephemerides->setEndDate("2100-01-01");
		ephemerides->generate();
		ephemerides->cancel();
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QCOMPARE(ephemerides->archivePath(), firstArchive);
		QCOMPARE(
		    QDir(exports.filePath("Ephemerides")).entryList(QDir::Dirs | QDir::Hidden | QDir::NoDotAndDotDot).size(),
		    1);
		ephemerides->setEndDate("2023-03-10");
		auto* table = ephemerides->table();
		auto* tableView = findItem(page, "ephemeris-table");
		QVERIFY(tableView);
		const auto columnOptions = table->columnOptions();
		for (const auto& option : columnOptions)
			table->setColumnEnabled(option.toMap().value("id").toInt(), false);
		QTRY_VERIFY_WITH_TIMEOUT(tableView->property("contentWidth").toDouble() >=
		                                 tableView->width() - tableView->property("rightMargin").toDouble() - 2 &&
		                             tableView->property("contentWidth").toDouble() <= tableView->width(),
		                         5000);
		for (const auto& option : columnOptions)
			table->setColumnEnabled(option.toMap().value("id").toInt(), true);
		QTRY_VERIFY_WITH_TIMEOUT(tableView->property("contentWidth").toDouble() > tableView->width(), 5000);
		for (const auto& option : columnOptions) {
			const auto values = option.toMap();
			table->setColumnEnabled(values.value("id").toInt(), values.value("enabled").toBool());
		}
		QCOMPARE(table->header(0), QString("Date UTC"));
		QCOMPARE(table->cell(0, 0), QString("2023-02-24 00:00"));
		QCOMPARE(table->cell(14, 0), QString("2023-03-10 00:00"));
		QVERIFY(table->cell(0, 5) != QString::fromUtf8("\u2014"));
		QVERIFY2(ephemerides->summary().startsWith("Geocentric (MPC 500)"), qPrintable(ephemerides->summary()));
		QVERIFY(ephemerides->summary().contains("two-body"));
		const auto originalRows = table->rows();
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("station500.json")), 1));
		QVERIFY(!ephemerides->busy());
		QCOMPARE(table->rows().front().longitude, originalRows.front().longitude);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QFile stationJson(exports.filePath("station500.json"));
		QVERIFY(stationJson.open(QIODevice::ReadOnly));
		const auto stationDocument = QJsonDocument::fromJson(stationJson.readAll()).object();
		QCOMPARE(stationDocument.value("center").toString(), QString("geocentric"));
		QCOMPARE(stationDocument.value("station").toString(), QString("500"));
		QCOMPARE(stationDocument.value("elongation_definition").toString(), QString("apparent"));
		QVERIFY(stationDocument.value("uncertainty_computed").toBool());
		QVERIFY(stationDocument.value("rows").toArray().first().toObject().value("range_3sigma_km").isDouble());
		QVERIFY(stationDocument.value("rows").toArray().first().toObject().value("range_rate_3sigma_km_s").isDouble());
		QVERIFY(!table->uncertaintySelected());
		click("ephemeris-settings-button");
		auto* ellipseOption = findItem(dialogContent, "ephemeris-column-ellipseMajor");
		QVERIFY(ellipseOption);
		QVERIFY(!ellipseOption->property("checked").toBool());
		QVERIFY(QMetaObject::invokeMethod(ellipseOption, "click"));
		QVERIFY(table->includesColumn(EphemerisTableModel::EllipseMajor));
		QVERIFY(table->includesColumn(EphemerisTableModel::EllipseMinor));
		QVERIFY(table->includesColumn(EphemerisTableModel::EllipseAngle));
		QVERIFY(QMetaObject::invokeMethod(ephemerisDialog, "close"));
		QVERIFY(!ephemerides->busy());
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QVERIFY(!ephemerides->uncertaintyUnavailable());
		QVERIFY(table->rows().front().uncertainty);
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("uncertainty.json")), 1));
		QFile uncertaintyJson(exports.filePath("uncertainty.json"));
		QVERIFY(uncertaintyJson.open(QIODevice::ReadOnly));
		const auto uncertaintyDocument = QJsonDocument::fromJson(uncertaintyJson.readAll()).object();
		QVERIFY(uncertaintyDocument.value("uncertainty_computed").toBool());
		const auto uncertaintyRow = uncertaintyDocument.value("rows").toArray().first().toObject();
		for (const auto* key : {"ra_cosdec_3sigma_arcsec", "dec_3sigma_arcsec", "ellipse_major_3sigma_arcsec",
		                        "ellipse_minor_3sigma_arcsec", "ellipse_angle_deg", "ellipse_area_3sigma_arcsec2",
		                        "sky_3sigma_arcsec", "range_3sigma_km", "range_rate_3sigma_km_s"})
			QVERIFY2(uncertaintyRow.value(QLatin1String(key)).isDouble(), key);
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("uncertainty.txt")), 0));
		QFile uncertaintyText(exports.filePath("uncertainty.txt"));
		QVERIFY(uncertaintyText.open(QIODevice::ReadOnly));
		const auto uncertaintyTextContent = uncertaintyText.readAll();
		for (const auto* label : {"RA cos Dec 3sigma", "Dec 3sigma", "Ellipse major 3sigma", "Ellipse minor 3sigma",
		                          "Ellipse angle", "Ellipse area 3sigma", "Sky 3sigma RSS", "Range 3sigma"})
			QVERIFY2(uncertaintyTextContent.contains(label), label);
		QCOMPARE(uncertaintyTextContent.count("Ellipse major 3sigma"), 1);
		table->setColumnEnabled(EphemerisTableModel::EllipseMajor, false);
		table->setColumnEnabled(EphemerisTableModel::Distance, false);
		table->setColumnEnabled(EphemerisTableModel::SkyMotion, true);
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("selected.txt")), 0));
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("selected.json")), 1));
		QFile selectedText(exports.filePath("selected.txt"));
		QVERIFY(selectedText.open(QIODevice::ReadOnly));
		const auto selectedTextContent = QString::fromUtf8(selectedText.readAll());
		QVERIFY(selectedTextContent.contains("Sky motion"));
		QVERIFY(!selectedTextContent.contains("Delta, AU"));
		QVERIFY(selectedTextContent.contains("Ellipse major 3sigma"));
		QVERIFY(selectedTextContent.contains("Range 3sigma"));
		QFile selectedJson(exports.filePath("selected.json"));
		QVERIFY(selectedJson.open(QIODevice::ReadOnly));
		const auto selectedRow =
		    QJsonDocument::fromJson(selectedJson.readAll()).object().value("rows").toArray().first().toObject();
		QVERIFY(selectedRow.contains("sky_motion_arcsec_h"));
		QVERIFY(!selectedRow.contains("delta_au"));
		QVERIFY(selectedRow.value("ellipse_major_3sigma_arcsec").isDouble());
		QVERIFY(selectedRow.value("range_3sigma_km").isDouble());
		QVERIFY(selectedRow.value("range_rate_3sigma_km_s").isDouble());
		table->setColumnEnabled(EphemerisTableModel::Distance, true);
		table->setColumnEnabled(EphemerisTableModel::SkyMotion, false);
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("spherical.txt")), 0));
		QFile spherical(exports.filePath("spherical.txt"));
		QVERIFY(spherical.open(QIODevice::ReadOnly));
		const auto sphericalLines = QString::fromUtf8(spherical.readAll()).split(QLatin1Char('\n'));
		const auto start = sphericalLines.indexOf(QStringLiteral("$$SOE"));
		QVERIFY(start >= 2);
		QVERIFY(sphericalLines[start - 2].contains("RA, h m s"));
		QVERIFY(sphericalLines[start - 2].contains("Delta, AU"));
		QVERIFY(sphericalLines.contains(QStringLiteral(
		    "# Elongation: apparent directions with light time, solar deflection and stellar aberration.")));
		QVERIFY(sphericalLines[start - 1].contains("---"));
		for (int row = 1; row <= ephemerides->epochCount(); ++row)
			QCOMPARE(sphericalLines[start + row].size(), sphericalLines[start - 2].size());
		const auto topocentric = table->rows();
		for (const auto& row : topocentric) {
			QCOMPARE(row.mjd, row.mjdUtc);
			QVERIFY(row.visualMagnitude.has_value());
			QVERIFY(row.distanceAu > 0.0 && row.heliocentricAu > 0.0);
		}

		click("ephemeris-center-option-1");
		QCOMPARE(ephemerides->center(), 1);
		QTRY_VERIFY(!findItem(page, "ephemeris-station")->isVisible());
		click("ephemeris-generate");
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy() && ephemerides->summary().startsWith("Geocentric"), 30000);
		QVERIFY(!ephemerides->archivePath().isEmpty());
		QVERIFY(ephemerides->archivePath() != firstArchive);
		QVERIFY(QFileInfo::exists(QDir(firstArchive).filePath("result.json")));
		for (std::size_t i = 0; i < topocentric.size(); ++i) {
			QCOMPARE(table->rows()[i].longitude, topocentric[i].longitude);
			QCOMPARE(table->rows()[i].latitude, topocentric[i].latitude);
		}

		click("ephemeris-layout-option-1");
		QCOMPARE(table->layout(), EphemerisTableModel::Layout::Cartesian);
		QCOMPARE(table->header(1), QString("X, AU"));
		QCOMPARE(ephemerides->epochCount(), 15);

		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("e.txt")), 0));
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("e.json")), 1));
		QVERIFY(!ephemerides->save(QUrl(QStringLiteral("https://example.com/e.json")), 1));
		QVERIFY(!ephemerides->error().isEmpty());
		QVERIFY(!ephemerides->save(QUrl::fromLocalFile(exports.filePath("missing/e.txt")), 0));
		QVERIFY(ephemerides->error().contains(QStringLiteral("e.txt")));
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("e.txt")), 0));
		QVERIFY(ephemerides->error().isEmpty());
		QFile json(exports.filePath("e.json"));
		QVERIFY(json.open(QIODevice::ReadOnly));
		const auto document = QJsonDocument::fromJson(json.readAll()).object();
		QCOMPARE(document.value("center").toString(), QString("geocentric"));
		QCOMPARE(document.value("time_scale").toString(), QString("UTC"));
		QCOMPARE(document.value("range_rate_definition").toString(), QString("radial_relative_velocity"));
		QVERIFY(document.value("solar_range_light_time").toBool());
		QCOMPARE(document.value("rows").toArray().size(), 15);
		QCOMPARE(document.value("layout").toString(), QString("cartesian"));
		QVERIFY(document.value("rows").toArray().first().toObject().value("x_au").isDouble());
		QVERIFY(!document.value("rows").toArray().first().toObject().contains("ra_deg"));
		QVERIFY(!document.value("rows").toArray().first().toObject().contains("delta_rate_km_s"));
		QFile text(exports.filePath("e.txt"));
		QVERIFY(text.open(QIODevice::ReadOnly));
		const QString textContent = QString::fromUtf8(text.readAll());
		QVERIFY(textContent.contains("Vx, AU/d"));
		QVERIFY(textContent.contains("light time applied"));
		QVERIFY(textContent.contains("not derivatives of light-time-corrected ranges"));
		QVERIFY(textContent.contains("$$SOE\n"));
		QVERIFY(textContent.contains("$$EOE\n"));
		const auto textLines = textContent.split(QLatin1Char('\n'));
		const auto textStart = textLines.indexOf(QStringLiteral("$$SOE"));
		QVERIFY(textStart >= 2);
		for (int row = 1; row <= ephemerides->epochCount(); ++row)
			QCOMPARE(textLines[textStart + row].size(), textLines[textStart - 2].size());
		text.close();
		table->setColumnEnabled(EphemerisTableModel::Distance, true);
		QCOMPARE(table->columns(), 8);
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("cartesian-extra.json")), 1));
		QFile cartesianExtra(exports.filePath("cartesian-extra.json"));
		QVERIFY(cartesianExtra.open(QIODevice::ReadOnly));
		const auto extraRow =
		    QJsonDocument::fromJson(cartesianExtra.readAll()).object().value("rows").toArray().first().toObject();
		QVERIFY(extraRow.contains("delta_au"));
		QVERIFY(extraRow.contains("x_au"));
		table->setColumnEnabled(EphemerisTableModel::Distance, false);
		table->setColumnEnabled(EphemerisTableModel::Elongation, true);
		ephemerides->setLightTime(false);
		ephemerides->generate();
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("geometric.txt")), 0));
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QFile geometric(exports.filePath("geometric.txt"));
		QVERIFY(geometric.open(QIODevice::ReadOnly));
		const auto geometricContent = QString::fromUtf8(geometric.readAll());
		QVERIFY(geometricContent.contains("geometric; no light time"));
		QVERIFY(!geometricContent.contains("Solar range includes Sun-to-target light time"));
		QVERIFY(geometricContent.contains("# Elongation: geometric Sun and target directions."));
		ephemerides->setLightTime(true);

		ephemerides->setScale(1);
		QCOMPARE(ephemerides->scaleName(), QString("TT"));
		click("ephemeris-generate");
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy() && table->header(0) == QString("Date TT"), 30000);
		QCOMPARE(table->rows().front().mjd, table->rows().front().mjdTt);
		QCOMPARE(table->rows().front().mjdTt, 59999.0);

		ephemerides->setEndDate("2023-02-01");
		QVERIFY(!ephemerides->formError().isEmpty());
		QVERIFY(!ephemerides->canGenerate());
		QTRY_VERIFY(!findItem(page, "ephemeris-generate")->isEnabled());
		ephemerides->setEndDate("2023-03-10");
		ephemerides->setStartTime("25:00");
		QVERIFY(!ephemerides->timeValid());
		QVERIFY(!ephemerides->canGenerate());
		ephemerides->setStartTime("00:00:00");
		ephemerides->setStepDays(0.0);
		QVERIFY(!ephemerides->canGenerate());
		ephemerides->setStepDays(1.0);

		ephemerides->setCenter(2);
		ephemerides->setStation("zzz");
		QCOMPARE(ephemerides->station(), QString("ZZZ"));
		QVERIFY(ephemerides->canGenerate());
		ephemerides->generate();
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY(!ephemerides->error().isEmpty());
		QCOMPARE(ephemerides->epochCount(), 15);
		QVERIFY(ephemerides->summary().startsWith("Geocentric"));

		ephemerides->setSolution(result, "Fixture", std::nullopt, nullptr, exports.path());
		QCOMPARE(ephemerides->epochCount(), 0);
		QVERIFY(ephemerides->error().isEmpty());
		QVERIFY(ephemerides->summary().isEmpty());
		QCOMPARE(ephemerides->startDate(), QString("2023-02-24"));
		ephemerides->setSolution(result, "C/2020 F3 #1", 18.0, nullptr, exports.path());
		auto withoutCovariance = result;
		withoutCovariance.fit.uncertainty = {};
		ephemerides->setSolution(withoutCovariance, "Fixture", 18.0, nullptr, exports.path());
		ephemerides->setCenter(1);
		ephemerides->generate();
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QVERIFY2(ephemerides->error().isEmpty(), qPrintable(ephemerides->error()));
		QVERIFY(ephemerides->save(QUrl::fromLocalFile(exports.filePath("no-covariance.json")), 1));
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QFile missingCovariance(exports.filePath("no-covariance.json"));
		QVERIFY(missingCovariance.open(QIODevice::ReadOnly));
		const auto missingDocument = QJsonDocument::fromJson(missingCovariance.readAll()).object();
		QVERIFY(!missingDocument.value("uncertainty_computed").toBool());
		QVERIFY(missingDocument.value("rows").toArray().first().toObject().value("range_3sigma_km").isNull());
		QVERIFY(missingDocument.value("rows").toArray().first().toObject().value("range_rate_3sigma_km_s").isNull());

		QFile invalidRoot(exports.filePath("not-a-directory"));
		QVERIFY(invalidRoot.open(QIODevice::WriteOnly));
		invalidRoot.close();
		ephemerides->setSolution(result, "Fixture", 18.0, nullptr, invalidRoot.fileName());
		ephemerides->generate();
		QTRY_VERIFY_WITH_TIMEOUT(!ephemerides->busy(), 30000);
		QCOMPARE(ephemerides->epochCount(), 15);
		QVERIFY(ephemerides->archivePath().isEmpty());
		QVERIFY(ephemerides->error().contains("could not be saved"));
	}

	void cleanupTestCase() { QSettings::setDefaultFormat(QSettings::NativeFormat); }
};
