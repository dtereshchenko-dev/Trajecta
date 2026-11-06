#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/FitConfig.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/ImpactRiskRmsObjective.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <vector>

namespace oc = orbital_core;

namespace {
	// Six C40 rows within one night (one vfcc17 batch of six) and two geocentric
	// rows days apart, JD-sorted.
	std::vector<ObsNode> observations()
	{
		std::vector<ObsNode> rows;
		const auto add = [&rows](double jd, const char* station, double ra) {
			ObsNode node;
			node.JD = jd;
			node.RA = ra;
			node.DEC = 0.1;
			node.obsCode = station;
			node.note2 = 'C';
			rows.push_back(node);
		};
		for (int i = 0; i < 6; ++i)
			add(2459165.5 + 0.01 * i, "C40", 1.0 + 0.001 * i);
		add(2459168.5, "500", 1.05);
		add(2459171.5, "500", 1.10);
		return rows;
	}

	oc::FitConfig vfcc17()
	{
		oc::FitConfig config;
		config.weightProfile = oc::AstrometricWeightProfile::Vfcc17;
		config.useWeights = true;
		return config;
	}
}  // namespace

TEST(ImpactRiskRmsObjectiveTest, SelectedRowsKeepTheFitRelaxationOfTheirBatch)
{
	const auto rows = observations();
	const Vec3D r(1.4, 0.4, 0.1);
	const Vec3D v(-0.004, 0.012, 0.001);
	const double epoch = rows.front().epochJdTdb();
	const auto orbit = make_orbit_from_states_vectors(StateVectors(r, v / Constants::GAUSS_K), epoch);
	ASSERT_TRUE(orbit);
	Eigen::VectorXd state(6);
	state << r.x, r.y, r.z, v.x, v.y, v.z;
	const NonGravModel noNonGrav;

	// The fit kept four of the six C40 rows and both geocentric rows.
	const std::vector<std::size_t> kept{0, 1, 4, 5, 6, 7};
	std::vector<ObsNode> keptRows;
	for (const auto index : kept)
		keptRows.push_back(rows[index]);
	const std::vector<ObsNode> keptC40{rows[0], rows[1], rows[4], rows[5]};

	const auto fitWeights = oc::hazard::makeImpactRiskRmsScaling(orbit, rows, kept, vfcc17(), noNonGrav);
	const auto keptOnly = oc::hazard::makeImpactRiskRmsScaling(orbit, keptRows, {}, vfcc17(), noNonGrav);
	const auto c40Only = oc::hazard::makeImpactRiskRmsScaling(orbit, keptC40, {}, vfcc17(), noNonGrav);
	ASSERT_TRUE(fitWeights.objective);
	ASSERT_TRUE(keptOnly.objective);
	ASSERT_TRUE(c40Only.objective);

	const double fit = fitWeights.objective(state);
	const double recounted = keptOnly.objective(state);
	const double c40 = c40Only.objective(state);
	ASSERT_TRUE(std::isfinite(fit) && std::isfinite(recounted) && std::isfinite(c40));

	// The fit counts the whole night, so the four kept C40 rows keep their sqrt(6/4)
	// relaxation and contribute 4/6 of their unrelaxed chi-square.  Counting only the
	// kept rows (a batch of four) would not relax them at all.
	EXPECT_GT(recounted - fit, 0.0);
	EXPECT_NEAR(fit, recounted - c40 * (1.0 - 4.0 / 6.0), 1.0e-9 * recounted);
}

TEST(ImpactRiskRmsObjectiveTest, EmptySelectionKeepsEveryRow)
{
	const auto rows = observations();
	const Vec3D r(1.4, 0.4, 0.1);
	const Vec3D v(-0.004, 0.012, 0.001);
	const double epoch = rows.front().epochJdTdb();
	const auto orbit = make_orbit_from_states_vectors(StateVectors(r, v / Constants::GAUSS_K), epoch);
	ASSERT_TRUE(orbit);
	Eigen::VectorXd state(6);
	state << r.x, r.y, r.z, v.x, v.y, v.z;
	const NonGravModel noNonGrav;

	std::vector<std::size_t> everyRow(rows.size());
	for (std::size_t i = 0; i < everyRow.size(); ++i)
		everyRow[i] = i;
	const auto explicitAll = oc::hazard::makeImpactRiskRmsScaling(orbit, rows, everyRow, vfcc17(), noNonGrav);
	const auto implicitAll = oc::hazard::makeImpactRiskRmsScaling(orbit, rows, {}, vfcc17(), noNonGrav);
	ASSERT_TRUE(explicitAll.objective);
	ASSERT_TRUE(implicitAll.objective);
	EXPECT_EQ(explicitAll.objective(state), implicitAll.objective(state));
}
