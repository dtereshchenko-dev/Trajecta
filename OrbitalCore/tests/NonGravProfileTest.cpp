#include <gtest/gtest.h>

#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/nongrav/NonGravProfile.h"

#include <vector>

TEST(NonGravProfileTest, NumFitParamsRoutesThroughModel)
{
	NonGravModel m;
	m.setEnabled(true);
	NonGravProfile profile(m);

	m.setModel(NgModel::None);
	EXPECT_EQ(profile.model().numFitParams(), 0);
	m.setModel(NgModel::Yarkovsky);
	EXPECT_EQ(profile.model().numFitParams(), 1);
	m.setModel(NgModel::InverseSquare);
	EXPECT_EQ(profile.model().numFitParams(), 3);
	m.setModel(NgModel::MarsdenH2O);
	EXPECT_EQ(profile.model().numFitParams(), 4);
	m.setModel(NgModel::InverseSquare);
	m.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 | NgParameterMask::DT);
	EXPECT_EQ(profile.model().numFitParams(), 4);
}

TEST(NonGravProfileTest, ExtendParamsAddsRightCount)
{
	NonGravModel m;
	m.setEnabled(true);
	NonGravProfile profile(m);

	{
		m.setModel(NgModel::InverseSquare);
		std::vector<double> p = {1, 2, 3, 4, 5, 6};
		std::vector<double> d = {1, 1, 1, 1, 1, 1};
		profile.extendParams(p, d);
		EXPECT_EQ(p.size(), 9u);
		EXPECT_EQ(d.size(), 9u);
	}
	{
		m.setModel(NgModel::MarsdenH2O);
		std::vector<double> p = {1, 2, 3, 4, 5, 6};
		std::vector<double> d = {1, 1, 1, 1, 1, 1};
		profile.extendParams(p, d);
		EXPECT_EQ(p.size(), 10u);
		EXPECT_EQ(d.size(), 10u);
		EXPECT_DOUBLE_EQ(p[9], 0.0);
	}
	{
		m.setModel(NgModel::Yarkovsky);
		std::vector<double> p = {1, 2, 3, 4, 5, 6};
		std::vector<double> d = {1, 1, 1, 1, 1, 1};
		profile.extendParams(p, d);
		EXPECT_EQ(p.size(), 7u);
		EXPECT_EQ(d.size(), 7u);
	}
}

TEST(NonGravProfileTest, ExtendParamsSeedsFromCurrentModel)
{
	NonGravModel m;
	m.setEnabled(true);
	m.setModel(NgModel::MarsdenH2O);
	m.setParams(1.5e-8, -2.5e-8, 0.7e-8);
	m.setDT(-30.0);

	NonGravProfile profile(m);

	std::vector<double> p, d;
	profile.extendParams(p, d);

	ASSERT_EQ(p.size(), 4u);
	EXPECT_DOUBLE_EQ(p[0], 1.5e-8);
	EXPECT_DOUBLE_EQ(p[1], -2.5e-8);
	EXPECT_DOUBLE_EQ(p[2], 0.7e-8);
	EXPECT_DOUBLE_EQ(p[3], -30.0);
}

TEST(NonGravProfileTest, UpdateFromParamsReadsDTForMarsden)
{
	NonGravModel m;
	m.setEnabled(true);
	m.setModel(NgModel::MarsdenH2O);
	NonGravProfile profile(m);

	std::vector<double> p = {0, 0, 0, 0, 0, 0, 1e-8, 2e-8, 3e-8, -40.0};
	profile.updateFromParams(p, 6);

	EXPECT_DOUBLE_EQ(m.A1(), 1e-8);
	EXPECT_DOUBLE_EQ(m.A2(), 2e-8);
	EXPECT_DOUBLE_EQ(m.A3(), 3e-8);
	EXPECT_DOUBLE_EQ(m.DT(), -40.0);
}

TEST(NonGravProfileTest, UpdateFromParamsPreservesFixedDTForInverseSquare)
{
	NonGravModel m;
	m.setEnabled(true);
	m.setModel(NgModel::InverseSquare);
	m.setDT(99.0);
	NonGravProfile profile(m);

	std::vector<double> p = {0, 0, 0, 0, 0, 0, 1e-8, 2e-8, 3e-8};
	profile.updateFromParams(p, 6);

	EXPECT_DOUBLE_EQ(m.A1(), 1e-8);
	EXPECT_DOUBLE_EQ(m.DT(), 99.0);
}

TEST(NonGravProfileTest, UpdateFromParamsForcesA1A3ZeroOnYarkovsky)
{
	NonGravModel m;
	m.setEnabled(true);
	m.setModel(NgModel::Yarkovsky);
	m.setParams(7.7e-8, 0.0, 5.5e-8);
	NonGravProfile profile(m);

	std::vector<double> p = {0, 0, 0, 0, 0, 0, -1.5e-14};
	profile.updateFromParams(p, 6);

	EXPECT_DOUBLE_EQ(m.A1(), 0.0);
	EXPECT_DOUBLE_EQ(m.A2(), -1.5e-14);
	EXPECT_DOUBLE_EQ(m.A3(), 0.0);
}

TEST(NonGravProfileTest, ParameterMaskPreservesFixedCoefficients)
{
	NonGravModel model;
	model.setEnabled(true);
	model.setModel(NgModel::MarsdenH2O);
	model.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2);
	model.setParams(1e-8, 2e-8, 3e-8);
	model.setDT(45.0);
	NonGravProfile profile(model);
	std::vector<double> params(6, 0.0), steps(6, 1.0);
	profile.extendParams(params, steps);
	ASSERT_EQ(params.size(), 8u);
	EXPECT_DOUBLE_EQ(params[6], 1e-8);
	EXPECT_DOUBLE_EQ(params[7], 2e-8);
	EXPECT_DOUBLE_EQ(steps[6], 1e-14);
	EXPECT_DOUBLE_EQ(steps[7], 1e-14);
	params[6] = 5e-8;
	params[7] = -4e-8;
	profile.updateFromParams(params, 6);
	EXPECT_DOUBLE_EQ(model.A1(), 5e-8);
	EXPECT_DOUBLE_EQ(model.A2(), -4e-8);
	EXPECT_DOUBLE_EQ(model.A3(), 3e-8);
	EXPECT_DOUBLE_EQ(model.DT(), 45.0);
}

TEST(NonGravProfileTest, DisabledModelDoesNotExtendTheFit)
{
	NonGravModel model;
	model.setModel(NgModel::MarsdenH2O);
	model.setEnabled(false);
	NonGravProfile profile(model);
	std::vector<double> params(6, 2.0), steps(6, 1.0);
	profile.extendParams(params, steps);
	EXPECT_EQ(params, std::vector<double>(6, 2.0));
	EXPECT_EQ(steps, std::vector<double>(6, 1.0));
}
