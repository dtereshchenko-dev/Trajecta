#ifndef TST_ORBITENGINE_H
#define TST_ORBITENGINE_H

#include <QTest>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QFile>
#include <QScopeGuard>
#include "controllers/OrbitEngine.h"
#include "controllers/OrbitWorker.h"
#include "settings/FitSettingsAggregator.h"
#include "settings/PerturbationsSettings.h"
#include <orbital_core/IEphemerisSource.h>
#include <orbital_core/DE440Source.h>
#include <cmath>
#include <limits>

class TestOrbitEngine : public QObject {
	Q_OBJECT

   private slots:

	void obsoleteIntegratorAccuracySettingsAreIgnoredWithoutCleanup()
	{
		QTemporaryDir dir;
		QVERIFY(dir.isValid());
		QSettings saved(dir.filePath("settings.ini"), QSettings::IniFormat);
		saved.setValue("integratorTolerance", 1e-3);
		saved.setValue("rk4Tolerance", 1e-2);
		saved.setValue("integrator/initialStepDays", 2.0);

		PerturbationsSettings settings;
		settings.load(saved);
		QCOMPARE(settings.integratorInitialStepDays(), 2.0);
		QCOMPARE(settings.metaObject()->indexOfProperty("integratorTolerance"), -1);
		settings.save(saved);
		QCOMPARE(saved.value("integratorTolerance").toDouble(), 1e-3);
		QCOMPARE(saved.value("rk4Tolerance").toDouble(), 1e-2);
		QCOMPARE(saved.value("integrator/initialStepDays").toDouble(), 2.0);
	}

	void fitControlsReachWorkerConfig()
	{
		OptimizationSettings settings;
		settings.forceModel.relativity = orbital_core::RelativityMode::SunOuter;
		settings.forceModel.earthJ3 = false;
		settings.forceModel.rotatingEarthPole = false;
		const auto forces = OrbitWorker::makeFitConfig(settings).nbody.forceModel;
		QVERIFY(forces.relativity == orbital_core::RelativityMode::SunOuter);
		QVERIFY(!forces.earthJ3);
		QVERIFY(!forces.rotatingEarthPole);
		settings.useWeights = true;
		QVERIFY(OrbitWorker::makeFitConfig(settings).useWeights);
		settings.useWeights = false;
		QVERIFY(!OrbitWorker::makeFitConfig(settings).useWeights);
	}

	void explicitSettingsAndInitialState()
	{
		FitSettingsAggregator settings;
		OrbitEngine engine(settings);
		QVERIFY(engine.errorMessage().isEmpty());
		QVERIFY(engine.progressStatus().isEmpty());
		QCOMPARE(settings.metaObject()->indexOfProperty("lm"), -1);
		QCOMPARE(settings.metaObject()->indexOfProperty("impactMonitor"), -1);
	}

	void solutionKeepsSettingsAtFitStart()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		DE440Source ephemeris;
		if (!ephemeris.isAvailable())
			QSKIP("DE440 is not available");
		auto* previousEphemeris = getDefaultEphemerisSource();
		setDefaultEphemerisSource(&ephemeris);
		const auto originalDir = getExecutableDir();
		const QByteArray originalXdg = qgetenv("XDG_DATA_HOME");
		const bool hadXdg = qEnvironmentVariableIsSet("XDG_DATA_HOME");
		const auto restore = qScopeGuard([&] {
			setDefaultEphemerisSource(previousEphemeris);
			getExecutableDir() = originalDir;
			if (hadXdg)
				qputenv("XDG_DATA_HOME", originalXdg);
			else
				qunsetenv("XDG_DATA_HOME");
		});
		getExecutableDir() = std::filesystem::path(directory.path().toStdWString());
		qputenv("XDG_DATA_HOME", directory.path().toUtf8());
		FitSettingsAggregator settings;
		settings.perturbations()->setProperty("earthJ3", false);
		OrbitEngine engine(settings);
		QSignalSpy completed(&engine, &OrbitEngine::solutionReady);
		QSignalSpy failed(&engine, &OrbitEngine::errorMessageChanged);
		engine.loadFile(QDir(QStringLiteral(ORBITAL_UI_SOURCE_DIR)).filePath("../data/examples/NEO/2020VV.xml"));
		settings.perturbations()->setProperty("earthJ3", true);
		QTRY_VERIFY_WITH_TIMEOUT(!completed.isEmpty() || !failed.isEmpty(), 120000);
		QVERIFY2(failed.isEmpty(), qPrintable(engine.errorMessage()));
		const QString folder = completed.first().at(1).toString();
		QSettings snapshot(QDir(folder).filePath("Trajecta.ini"), QSettings::IniFormat);
		QVERIFY(QFile::exists(snapshot.fileName()));
		QVERIFY(snapshot.contains("forces/earthJ3"));
		QVERIFY(!snapshot.value("forces/earthJ3").toBool());
		QVERIFY(QSettings().value("forces/earthJ3").toBool());
	}

	void loadFailureAndCancellation()
	{
		FitSettingsAggregator settings;
		OrbitEngine engine(settings);
		QSignalSpy failed(&engine, &OrbitEngine::errorMessageChanged);
		QSignalSpy completed(&engine, &OrbitEngine::solutionReady);
		engine.loadFile("nonexistent-observations.xml");
		QTRY_COMPARE_WITH_TIMEOUT(failed.count(), 1, 10000);
		QVERIFY(!engine.errorMessage().isEmpty());
		QCOMPARE(completed.count(), 0);

		engine.loadFile("nonexistent-observations.xml");
		engine.cancel();
		QTRY_COMPARE_WITH_TIMEOUT(failed.count(), 2, 10000);
		QCOMPARE(completed.count(), 0);
	}
};

#endif  // TST_ORBITENGINE_H
