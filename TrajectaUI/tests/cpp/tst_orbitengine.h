#ifndef TST_ORBITENGINE_H
#define TST_ORBITENGINE_H

#include <QTest>
#include <QSignalSpy>
#include <QTemporaryDir>
#include "controllers/OrbitEngine.h"
#include "controllers/OrbitWorker.h"
#include "settings/FitSettingsAggregator.h"
#include "settings/PerturbationsSettings.h"
#include <orbital_core/IEphemerisSource.h>
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
		saved.setValue("integratorInitialStepDays", 2.0);

		PerturbationsSettings settings;
		settings.load(saved);
		QCOMPARE(settings.integratorInitialStepDays(), 2.0);
		QCOMPARE(settings.metaObject()->indexOfProperty("integratorTolerance"), -1);
		settings.save(saved);
		QCOMPARE(saved.value("integratorTolerance").toDouble(), 1e-3);
		QCOMPARE(saved.value("rk4Tolerance").toDouble(), 1e-2);
		QCOMPARE(saved.value("integratorInitialStepDays").toDouble(), 2.0);
	}

	void fitControlsReachWorkerConfig()
	{
		OptimizationSettings settings;
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
