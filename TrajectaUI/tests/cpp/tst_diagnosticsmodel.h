#pragma once

#include <QTest>
#include <initializer_list>
#include "models/DiagnosticsModel.h"
#include "features/solution/SolutionViewModel.h"
#include "features/observations/ObservationReviewModel.h"

class TestDiagnosticsModel : public QObject {
	Q_OBJECT

	static Residuals residuals(std::initializer_list<double> values)
	{
		Residuals result;
		for (const double value : values) {
			Residual row{};
			row.dRA_arcsec = value;
			row.total_arcsec = value;
			result.add(row);
		}
		return result;
	}

   private slots:

	void emptyRms()
	{
		DiagnosticsModel model;
		model.setFromFit({}, 0.0, 0.0);
		QCOMPARE(model.rms(), 0.0);
	}

	void rmsRetainsAllResiduals()
	{
		const auto input = residuals({1.0, 1.0, 100.0});
		DiagnosticsModel model;
		model.setFromFit(input, 0.0, 0.0);
		QCOMPARE(model.rms(), input.rms());
		QCOMPARE(model.fillResidualCharts(nullptr, nullptr).value("count").toInt(), 3);
	}

	void clearResetsRms()
	{
		DiagnosticsModel model;
		model.setFromFit(residuals({2.0}), 0.0, 0.0);
		model.clear();
		QCOMPARE(model.rms(), 0.0);
	}

	void reviewBulkSelectionRespectsFilterAndReplacedObservations()
	{
		using namespace orbital_core;
		ObservationSet observations;
		observations.all.resize(4);
		for (std::size_t i = 0; i < observations.all.size(); ++i) {
			auto& o = observations.all[i];
			o.id = std::to_string(i);
			o.JD = 2451545.0;
			o.obsCode = "500";
		}
		observations.all[1].note2 = 'X';
		observations.coincidentObservations.push_back({{"0", "1"}, {ObservationIssue::CoincidentPosition}, 0.0});
		observations.visibilityDiagnostics.push_back({"2", {ObservationIssue::BelowHorizon}, -20.0, -30.0});
		ObservationReviewModel model;
		model.setObservations(observations);
		QCOMPARE(model.totalCount(), 3);
		model.setFilter(model.positionFilter());
		model.setAllShownIncluded(false);
		QCOMPARE(model.excludedCount(), 1);
		QCOMPARE(model.overrides().size(), std::size_t{1});
		QCOMPARE(QString::fromStdString(model.overrides()[0].id), QStringLiteral("0"));
		model.setFilter(model.belowHorizonFilter());
		model.setAllShownIncluded(false);
		QCOMPARE(model.excludedCount(), 2);
		model.setAllShownIncluded(true);
		QCOMPARE(model.excludedCount(), 1);
		model.setFilter(0);
		model.setAllShownIncluded(true);
		QCOMPARE(model.excludedCount(), 0);
	}

	void solutionUsesFittedRowsWithoutDiscardingInput()
	{
		orbital_core::FitPipelineResult result;
		result.observations.all.resize(3);
		result.observations.all[0].JD = 2459000.0;
		result.observations.all[1].JD = 2459010.0;
		result.observations.all[2].JD = 2459012.0;
		result.fit.keptObsIndices = {1, 2};
		result.fit.residuals = residuals({1.0, 2.0});
		result.fit.allResiduals = residuals({100.0, 1.0, 2.0});
		result.fit.orbit =
		    make_orbit_from_states_vectors(StateVectors(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.9, 0.1)), 2459010.0);
		SolutionViewModel solution;
		QVERIFY(solution.publish(result, "fixture"));
		QCOMPARE(solution.diagnostics()->arcLength(), 2.0);
		QCOMPARE(solution.diagnostics()->rms(), result.fit.residuals.rms());
		QCOMPARE(result.observations.all.size(), std::size_t{3});
		QCOMPARE(result.fit.keptObsIndices, (std::vector<std::size_t>{1, 2}));

		QVERIFY(!solution.publish({}, "failed"));
		QCOMPARE(solution.diagnostics()->arcLength(), 2.0);
		QCOMPARE(solution.diagnostics()->rms(), result.fit.residuals.rms());
	}
};
