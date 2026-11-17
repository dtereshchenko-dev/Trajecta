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
		stored.beginGroup("non-gravitational");

		NonGravSettings settings;
		QVERIFY(settings.setProperty("enabled", true));
		QVERIFY(settings.setProperty("model", 3));
		QVERIFY(settings.setProperty("solveMask", 15));
		settings.save(stored);
		QCOMPARE(stored.allKeys(), QStringList({"enabled", "model", "solveMask"}));
		stored.sync();
		QCOMPARE(stored.status(), QSettings::NoError);

		QSettings persisted(directory.filePath("settings.ini"), QSettings::IniFormat);
		persisted.beginGroup("non-gravitational");
		NonGravSettings reloaded;
		reloaded.load(persisted);
		QVERIFY(reloaded.enabled());
		QCOMPARE(reloaded.model(), 3);
		QCOMPARE(reloaded.solveMask(), 15);
		stored.setValue("solveMask", 0);
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

	void loadingAndSavingSettingsRemovesObsoleteKeys()
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
		QVERIFY(!stored.childGroups().contains("optimizer"));
	}

	void migratesSettingsToTheirSections()
	{
		QSettings stored;
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		stored.clear();
		stored.setValue("optimizer/pertEnabled", false);
		stored.setValue("optimizer/integratorInitialStepDays", 2.5);
		stored.setValue("optimizer/planetMask", 123u);
		stored.setValue("optimizer/useSb441", true);
		stored.setValue("optimizer/nonGravEnabled", true);
		stored.setValue("optimizer/nonGravModel", 3);
		stored.setValue("optimizer/nonGravSolveMask", 2);
		stored.setValue("optimizer/useWeights", false);
		stored.setValue("optimizer/sigmaClipThreshold", 4.0);
		stored.setValue("optimizer/photometry/asteroidG", 0.24);
		stored.setValue("optimizer/unusedSetting", 42);
		stored.sync();

		FitSettingsAggregator settings;
		settings.load();
		QVERIFY(!settings.perturbations()->enabled());
		QCOMPARE(settings.perturbations()->integratorInitialStepDays(), 2.5);
		QCOMPARE(settings.perturbations()->planetMask(), 123u);
		QVERIFY(settings.perturbations()->useSb441());
		QVERIFY(settings.nonGrav()->enabled());
		QCOMPARE(settings.nonGrav()->model(), 3);
		QCOMPARE(settings.nonGrav()->solveMask(), 2);
		QVERIFY(!settings.weighting()->useWeights());
		stored.sync();
		QCOMPARE(stored.value("fit/sigmaClipThreshold").toDouble(), 4.0);
		QCOMPARE(stored.value("photometry/asteroidG").toDouble(), 0.24);
		QVERIFY(!stored.childGroups().contains("optimizer"));

		settings.save();
		stored.sync();
		QVERIFY(!stored.value("integrator/nBodyPropagation").toBool());
		QCOMPARE(stored.value("perturbations/planetMask").toUInt(), 123u);
		QCOMPARE(stored.value("non-gravitational/model").toInt(), 3);
	}

	void keepsNewSettingsWhenOldValuesConflict()
	{
		QSettings stored;
		const auto clear = qScopeGuard([] { QSettings().clear(); });
		stored.clear();
		stored.setValue("fit/useWeights", false);
		stored.setValue("optimizer/useWeights", true);
		stored.sync();
		FitSettingsAggregator settings;
		settings.load();
		stored.sync();
		QVERIFY(!settings.weighting()->useWeights());
		QVERIFY(!stored.childGroups().contains("optimizer"));
	}
};
