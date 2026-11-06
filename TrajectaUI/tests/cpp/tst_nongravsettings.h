#pragma once

#include <QSettings>
#include <QTemporaryDir>
#include <QTest>
#include <QScopeGuard>

#include "FitSettingsAggregator.h"
#include "NonGravSettings.h"

class TestNonGravSettings : public QObject {
	Q_OBJECT

   private slots:

	void modelAndSolveMaskRoundTrip()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings stored(directory.filePath("settings.ini"), QSettings::IniFormat);
		stored.beginGroup("optimizer");

		NonGravSettings settings;
		QVERIFY(settings.setProperty("enabled", true));
		QVERIFY(settings.setProperty("model", 3));
		QVERIFY(settings.setProperty("solveMask", 15));
		settings.save(stored);
		QCOMPARE(stored.allKeys(), QStringList({"nonGravEnabled", "nonGravModel", "nonGravSolveMask"}));
		stored.sync();
		QCOMPARE(stored.status(), QSettings::NoError);

		QSettings persisted(directory.filePath("settings.ini"), QSettings::IniFormat);
		persisted.beginGroup("optimizer");
		NonGravSettings reloaded;
		reloaded.load(persisted);
		QVERIFY(reloaded.enabled());
		QCOMPARE(reloaded.model(), 3);
		QCOMPARE(reloaded.solveMask(), 15);
		stored.setValue("nonGravSolveMask", 0);
		reloaded.load(stored);
		QCOMPARE(reloaded.model(), 3);
		QCOMPARE(reloaded.solveMask(), 0);
	}

	void obsoleteValuesAreIgnoredWithoutModifyingStorage()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		QSettings stored(directory.filePath("settings.ini"), QSettings::IniFormat);
		for (const auto* key : {"a1Prior", "a2Prior", "a3Prior", "dtPrior", "a1Std", "a2Std", "a3Std", "dtStd",
		                        "nonGravRegStrength", "nonGravUsePriorAsInitial"})
			stored.setValue(key, 42.0);
		const auto keys = stored.allKeys();

		NonGravSettings settings;
		settings.load(stored);
		QVERIFY(!settings.enabled());
		QCOMPARE(settings.model(), 1);
		QCOMPARE(settings.solveMask(), 0);
		QCOMPARE(stored.allKeys(), keys);
		settings.save(stored);
		for (const auto& key : keys)
			QCOMPARE(stored.value(key).toDouble(), 42.0);
		for (const auto* property : {"a1Prior", "a2Prior", "a3Prior", "dtPrior", "a1Std", "a2Std", "a3Std", "dtStd",
		                             "regStrength", "usePriorAsInitial"})
			QCOMPARE(settings.metaObject()->indexOfProperty(property), -1);
	}

	void loadingAndSavingSettingsDoesNotCleanObsoleteKeys()
	{
		QSettings stored;
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		stored.clear();
		const QStringList keys{
		    "optimizer/lm/lambdaInit",     "optimizer/perturbations/step",  "optimizer/nonGrav/a2Prior",
		    "optimizer/weighting/profile", "optimizer/sigmaClip/threshold", "optimizer/impactMonitor/years",
		    "optimizer/outputEpochJD",     "optimizer/maxIterations",       "optimizer/lambdaInit",
		    "optimizer/lambdaUp",          "optimizer/lambdaDown",          "optimizer/integratorTolerance",
		    "optimizer/rk4StepSize",       "optimizer/rk4Tolerance",        "optimizer/overobsTimeSpanDay",
		    "optimizer/overobsCeiling",    "optimizer/timingSigmaSeconds"};
		for (const auto& key : keys)
			stored.setValue(key, 42.0);
		stored.sync();

		FitSettingsAggregator settings;
		settings.load();
		QCOMPARE(settings.perturbations()->integratorInitialStepDays(), 10.0);
		QVERIFY(!settings.weighting()->timingEnabled());
		settings.save();
		stored.sync();
		for (const auto& key : keys)
			QCOMPARE(stored.value(key).toDouble(), 42.0);
	}
};
