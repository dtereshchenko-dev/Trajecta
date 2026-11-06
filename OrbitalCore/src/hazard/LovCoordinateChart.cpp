#include "hazard/LovCoordinateChart.h"
#include "orbital_core/hazard/EquinoctialCoordinates.h"
#include "orbital_core/CometaryCoordinates.h"
#include <optional>

#include "orbital_core/Constants.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		using VecX = Eigen::VectorXd;
		using MatX = Eigen::MatrixXd;

		constexpr int kStateSize = 6;

		[[nodiscard]] VecX physicalParameters(const Orbit& orbit, const NonGravModel& ng)
		{
			const int nExtra = ng.numFitParams();
			VecX x(kStateSize + nExtra);
			const Vec3D r = orbit.getR();
			const Vec3D v = orbit.getV() * Constants::GAUSS_K;
			x(0) = r.x;
			x(1) = r.y;
			x(2) = r.z;
			x(3) = v.x;
			x(4) = v.y;
			x(5) = v.z;

			for (int k = 0; k < nExtra; ++k)
				x(6 + k) = ng.fitParameterValue(k);
			return x;
		}

		[[nodiscard]] std::optional<VecX> keplerianFromPhysical(const VecX& physical, double epochTdb,
		                                                        const NonGravModel& ng)
		{
			if (physical.size() < kStateSize)
				return std::nullopt;

			const Vec3D r(physical(0), physical(1), physical(2));
			const Vec3D vPhysical(physical(3), physical(4), physical(5));
			auto orbit =
			    make_orbit_from_states_vectors(StateVectors(r, vPhysical * (1.0 / Constants::GAUSS_K)), epochTdb);
			if (!orbit || orbit->getType() != OrbitType::Elliptical)
				return std::nullopt;

			const int nExtra = ng.numFitParams();
			VecX out(kStateSize + nExtra);
			out(0) = orbit->getSemiMajorAxis();
			out(1) = orbit->getExcentricity();
			out(2) = orbit->getInclination();
			out(3) = Helpers::wrapTwoPi(orbit->getLongitudeOfAscendingNode());
			out(4) = Helpers::wrapTwoPi(orbit->getArgumentOfPericentre());
			out(5) = Helpers::wrapTwoPi(orbit->getMeanAnomaly());
			for (int j = 0; j < nExtra; ++j)
				out(kStateSize + j) = physical(kStateSize + j);
			return out;
		}

		[[nodiscard]] std::optional<VecX> equinoctialFromPhysical(const VecX& physical, double, const NonGravModel& ng)
		{
			if (physical.size() != kStateSize + ng.numFitParams() || !physical.allFinite())
				return std::nullopt;
			const auto elements = equinoctial::fromCartesian(physical.head<6>());
			if (!elements)
				return std::nullopt;
			VecX out = physical;
			out.head<6>() = *elements;
			return out;
		}

		[[nodiscard]] std::optional<VecX> physicalFromKeplerian(const VecX& keplerian, double, const NonGravModel& ng)
		{
			if (keplerian.size() < kStateSize)
				return std::nullopt;

			const double a = keplerian(0);
			const double e = keplerian(1);
			const double incl = keplerian(2);
			const double node = keplerian(3);
			const double peri = keplerian(4);
			const double mean = keplerian(5);
			if (!(a > 0.0) || !(e >= 0.0) || !(e < 1.0) || !std::isfinite(incl) || !std::isfinite(node) ||
			    !std::isfinite(peri) || !std::isfinite(mean)) {
				return std::nullopt;
			}

			const double E = Calculations::kepler_solve_e(Helpers::wrapTwoPi(mean), e);
			const double cosE = std::cos(E);
			const double sinE = std::sin(E);
			const double oneMinusECosE = 1.0 - e * cosE;
			if (!(oneMinusECosE > 0.0))
				return std::nullopt;

			const double beta = std::sqrt(std::max(0.0, 1.0 - e * e));
			const double xOrb = a * (cosE - e);
			const double yOrb = a * beta * sinE;
			const double speedScale = 1.0 / (std::sqrt(a) * oneMinusECosE);
			const double vxOrb = -sinE * speedScale;
			const double vyOrb = beta * cosE * speedScale;

			const double cO = std::cos(node);
			const double sO = std::sin(node);
			const double ci = std::cos(incl);
			const double si = std::sin(incl);
			const double cw = std::cos(peri);
			const double sw = std::sin(peri);

			const std::array<double, 9> rot{cO * cw - sO * sw * ci,
			                                -cO * sw - sO * cw * ci,
			                                sO * si,
			                                sO * cw + cO * sw * ci,
			                                -sO * sw + cO * cw * ci,
			                                -cO * si,
			                                sw * si,
			                                cw * si,
			                                ci};

			const int nExtra = ng.numFitParams();
			VecX out(kStateSize + nExtra);
			out(0) = rot[0] * xOrb + rot[1] * yOrb;
			out(1) = rot[3] * xOrb + rot[4] * yOrb;
			out(2) = rot[6] * xOrb + rot[7] * yOrb;
			out(3) = Constants::GAUSS_K * (rot[0] * vxOrb + rot[1] * vyOrb);
			out(4) = Constants::GAUSS_K * (rot[3] * vxOrb + rot[4] * vyOrb);
			out(5) = Constants::GAUSS_K * (rot[6] * vxOrb + rot[7] * vyOrb);

			for (int j = 0; j < nExtra; ++j)
				out(kStateSize + j) = keplerian(kStateSize + j);
			return out;
		}

		[[nodiscard]] std::optional<VecX> physicalFromEquinoctial(const VecX& coordinates, double,
		                                                          const NonGravModel& ng)
		{
			if (coordinates.size() != kStateSize + ng.numFitParams() || !coordinates.allFinite())
				return std::nullopt;
			const auto state = equinoctial::toCartesian(coordinates.head<6>());
			if (!state)
				return std::nullopt;
			VecX out = coordinates;
			out.head<6>() = *state;
			return out;
		}

		[[nodiscard]] double finiteDifferenceStep(double x, int index) noexcept
		{
			const double absX = std::abs(x);
			const double floor = index < 3 ? 1.0e-8 : index < 6 ? 1.0e-10 : 1.0e-14;
			return std::max(floor, absX * 1.0e-7);
		}

		[[nodiscard]] double chartFiniteDifferenceStep(double x, int index) noexcept
		{
			if (index >= kStateSize)
				return 1.0;
			const double floor = index == 0 ? 1.0e-8 : 1.0e-9;
			return std::max(floor, std::abs(x) * 1.0e-7);
		}

		[[nodiscard]] std::optional<VecX> coordinatesFromPhysical(LovCoordinateSystem coordinateSystem,
		                                                          const VecX& physical, double epochTdb,
		                                                          const NonGravModel& ng)
		{
			if (coordinateSystem == LovCoordinateSystem::Keplerian)
				return keplerianFromPhysical(physical, epochTdb, ng);
			if (coordinateSystem == LovCoordinateSystem::Cometary) {
				if (physical.size() != kStateSize + ng.numFitParams() || !physical.allFinite())
					return std::nullopt;
				const auto elements = cometary::fromCartesian(physical.head<6>());
				if (!elements)
					return std::nullopt;
				VecX out = physical;
				out.head<6>() = *elements;
				return out;
			}
			return equinoctialFromPhysical(physical, epochTdb, ng);
		}

		[[nodiscard]] std::optional<VecX> physicalFromCoordinates(LovCoordinateSystem coordinateSystem,
		                                                          const VecX& coordinates, double epochTdb,
		                                                          const NonGravModel& ng)
		{
			if (coordinateSystem == LovCoordinateSystem::Keplerian)
				return physicalFromKeplerian(coordinates, epochTdb, ng);
			if (coordinateSystem == LovCoordinateSystem::Cometary) {
				if (coordinates.size() != kStateSize + ng.numFitParams() || !coordinates.allFinite())
					return std::nullopt;
				const auto state = cometary::toCartesian(coordinates.head<6>());
				if (!state)
					return std::nullopt;
				VecX out = coordinates;
				out.head<6>() = *state;
				return out;
			}
			return physicalFromEquinoctial(coordinates, epochTdb, ng);
		}

		[[nodiscard]] VecX coordinateDifference(const VecX& plus, const VecX& minus,
		                                        LovCoordinateSystem coordinateSystem)
		{
			VecX diff = plus - minus;
			if (diff.size() < kStateSize)
				return diff;

			if (coordinateSystem != LovCoordinateSystem::Equinoctial) {
				diff(3) = Helpers::wrapPi(plus(3) - minus(3));
				diff(4) = Helpers::wrapPi(plus(4) - minus(4));
				diff(5) = Helpers::wrapPi(plus(5) - minus(5));
			} else {
				diff(5) = Helpers::wrapPi(plus(5) - minus(5));
			}
			return diff;
		}

		[[nodiscard]] std::optional<MatX> dCoordinatesDPhysical(const VecX& physical, double epochTdb,
		                                                        const NonGravModel& ng,
		                                                        LovCoordinateSystem coordinateSystem)
		{
			const auto base = coordinatesFromPhysical(coordinateSystem, physical, epochTdb, ng);
			if (!base)
				return std::nullopt;

			const int n = static_cast<int>(physical.size());
			if (coordinateSystem != LovCoordinateSystem::Keplerian) {
				const auto orbitalJacobian = coordinateSystem == LovCoordinateSystem::Cometary
				                                 ? cometary::cartesianToElementsJacobian(base->head<6>())
				                                 : equinoctial::cartesianToElementsJacobian(base->head<6>());
				if (!orbitalJacobian)
					return std::nullopt;
				MatX jac = MatX::Identity(n, n);
				jac.topLeftCorner<6, 6>() = *orbitalJacobian;
				return jac;
			}
			MatX jac = MatX::Zero(n, n);
			for (int j = 0; j < n; ++j) {
				const double h = finiteDifferenceStep(physical(j), j);
				VecX xp = physical;
				VecX xm = physical;
				xp(j) += h;
				xm(j) -= h;

				const auto yp = coordinatesFromPhysical(coordinateSystem, xp, epochTdb, ng);
				const auto ym = coordinatesFromPhysical(coordinateSystem, xm, epochTdb, ng);
				if (!yp || !ym)
					return std::nullopt;

				jac.col(j) = coordinateDifference(*yp, *ym, coordinateSystem) * (0.5 / h);
			}
			return jac;
		}

		[[nodiscard]] std::optional<MatX> dPhysicalDCoordinates(const VecX& coordinates, double epochTdb,
		                                                        const NonGravModel& ng,
		                                                        LovCoordinateSystem coordinateSystem)
		{
			const auto base = physicalFromCoordinates(coordinateSystem, coordinates, epochTdb, ng);
			if (!base)
				return std::nullopt;

			const int n = static_cast<int>(coordinates.size());
			if (coordinateSystem != LovCoordinateSystem::Keplerian) {
				equinoctial::Matrix6 orbitalJacobian;
				const auto state = coordinateSystem == LovCoordinateSystem::Cometary
				                       ? cometary::toCartesian(coordinates.head<6>(), &orbitalJacobian)
				                       : equinoctial::toCartesian(coordinates.head<6>(), &orbitalJacobian);
				if (!state)
					return std::nullopt;
				MatX jac = MatX::Identity(n, n);
				jac.topLeftCorner<6, 6>() = orbitalJacobian;
				return jac;
			}
			MatX jac = MatX::Zero(n, n);
			for (int j = 0; j < n; ++j) {
				if (j >= kStateSize) {
					jac(j, j) = 1.0;
					continue;
				}

				const double h = chartFiniteDifferenceStep(coordinates(j), j);
				VecX yp = coordinates;
				VecX ym = coordinates;
				yp(j) += h;
				ym(j) -= h;

				const auto xp = physicalFromCoordinates(coordinateSystem, yp, epochTdb, ng);
				const auto xm = physicalFromCoordinates(coordinateSystem, ym, epochTdb, ng);
				if (!xp || !xm)
					return std::nullopt;

				jac.col(j) = (*xp - *xm) * (0.5 / h);
			}
			return jac;
		}

		[[nodiscard]] LovCoordinateChart makeLovChart(const Orbit& orbit,
		                                              const FitParameterCovariance& physicalCovariance,
		                                              const NonGravModel& nonGrav, LovCoordinateSystem coordinateSystem)
		{
			LovCoordinateChart chart;
			chart.nonGrav = nonGrav;
			chart.epochTdb = orbit.getEpoch();
			chart.coordinateSystem = coordinateSystem;

			const VecX physical = physicalParameters(orbit, nonGrav);
			if (physicalCovariance.rows() != physical.size() || physicalCovariance.cols() != physical.size() ||
			    !physicalCovariance.array().isFinite().all()) {
				return chart;
			}

			const auto coordinates = coordinatesFromPhysical(coordinateSystem, physical, orbit.getEpoch(), nonGrav);
			const auto jac = dCoordinatesDPhysical(physical, orbit.getEpoch(), nonGrav, coordinateSystem);
			if (!coordinates || !jac)
				return chart;

			chart.nominalCoordinates = *coordinates;
			chart.coordinateCovariance = (*jac) * physicalCovariance * jac->transpose();
			chart.coordinateCovariance =
			    0.5 * (chart.coordinateCovariance + chart.coordinateCovariance.transpose()).eval();
			chart.valid =
			    chart.coordinateCovariance.array().isFinite().all() && !chart.coordinateCovariance.isZero(0.0);
			return chart;
		}
	}  // namespace

	Eigen::VectorXd lovCoordinateUnits(LovCoordinateSystem coordinateSystem, const Eigen::VectorXd& coordinates,
	                                   bool scaling)
	{
		if (coordinates.size() < 6)
			return {};

		Eigen::VectorXd units = Eigen::VectorXd::Ones(coordinates.size());
		if (!scaling)
			return units;

		if (coordinateSystem == LovCoordinateSystem::Keplerian) {
			units(0) = 1.0e-3;
			units(2) = Constants::PI;
			units.segment<3>(3).setConstant(Constants::TWOPI);
		} else if (coordinateSystem == LovCoordinateSystem::Cometary) {
			units(0) = coordinates(0);
			units(2) = Constants::PI;
			units.segment<3>(3).setConstant(Constants::TWOPI);
		} else {
			units(0) = coordinates(0);
			units(5) = Constants::TWOPI;
		}
		if (!units.array().isFinite().all() || (units.array() == 0.0).any()) {
			return {};
		}
		return units;
	}

	LovCoordinateChart makeEquinoctialLovChart(const Orbit& orbit, const FitParameterCovariance& physicalCovariance,
	                                           const NonGravModel& nonGrav)
	{
		return makeLovChart(orbit, physicalCovariance, nonGrav, LovCoordinateSystem::Equinoctial);
	}

	LovCoordinateChart makeKeplerianLovChart(const Orbit& orbit, const FitParameterCovariance& physicalCovariance,
	                                         const NonGravModel& nonGrav)
	{
		return makeLovChart(orbit, physicalCovariance, nonGrav, LovCoordinateSystem::Keplerian);
	}

	std::optional<Eigen::VectorXd> keplerianCoordinates(const Orbit& orbit, const NonGravModel& nonGrav)
	{
		const VecX physical = physicalParameters(orbit, nonGrav);
		return keplerianFromPhysical(physical, orbit.getEpoch(), nonGrav);
	}

	LovCoordinateChart makeCometaryLovChart(const Orbit& orbit, const FitParameterCovariance& physicalCovariance,
	                                        const NonGravModel& nonGrav)
	{
		return makeLovChart(orbit, physicalCovariance, nonGrav, LovCoordinateSystem::Cometary);
	}

	LovCoordinateChart makeAutomaticLovChart(const Orbit& orbit, const FitParameterCovariance& physicalCovariance,
	                                         const NonGravModel& nonGrav, LovCoordinateSystem preferred)
	{
		const auto system =
		    cometary::usesCometaryElements(orbit.getExcentricity()) ? LovCoordinateSystem::Cometary : preferred;
		return makeLovChart(orbit, physicalCovariance, nonGrav, system);
	}

	std::optional<Eigen::VectorXd> cometaryCoordinates(const Orbit& orbit, const NonGravModel& nonGrav)
	{
		return coordinatesFromPhysical(LovCoordinateSystem::Cometary, physicalParameters(orbit, nonGrav),
		                               orbit.getEpoch(), nonGrav);
	}

	std::optional<Eigen::VectorXd> equinoctialCoordinates(const Orbit& orbit, const NonGravModel& nonGrav)
	{
		return equinoctialFromPhysical(physicalParameters(orbit, nonGrav), orbit.getEpoch(), nonGrav);
	}

	std::optional<Eigen::VectorXd> physicalParametersFromChart(const LovCoordinateChart& chart,
	                                                           const Eigen::VectorXd& coordinates)
	{
		if (!chart.valid || coordinates.size() != chart.nominalCoordinates.size())
			return std::nullopt;

		return physicalFromCoordinates(chart.coordinateSystem, coordinates, chart.epochTdb, chart.nonGrav);
	}

	bool physicalParametersFromChart(const LovCoordinateChart& chart, const Eigen::VectorXd& coordinates,
	                                 Eigen::Ref<Eigen::VectorXd> output)
	{
		if (output.size() != coordinates.size())
			return false;

		const auto physical = physicalParametersFromChart(chart, coordinates);
		if (!physical || physical->size() != output.size())
			return false;

		output = *physical;
		return output.array().isFinite().all();
	}

	std::optional<Eigen::MatrixXd> physicalJacobianFromChart(const LovCoordinateChart& chart,
	                                                         const Eigen::VectorXd& coordinates)
	{
		if (!chart.valid || coordinates.size() != chart.nominalCoordinates.size()) {
			return std::nullopt;
		}

		return dPhysicalDCoordinates(coordinates, chart.epochTdb, chart.nonGrav, chart.coordinateSystem);
	}

	bool physicalJacobianFromChart(const LovCoordinateChart& chart, const Eigen::VectorXd& coordinates,
	                               Eigen::Ref<Eigen::MatrixXd> output)
	{
		if (output.rows() != coordinates.size() || output.cols() != coordinates.size()) {
			return false;
		}

		const auto jacobian = physicalJacobianFromChart(chart, coordinates);
		if (!jacobian || jacobian->rows() != output.rows() || jacobian->cols() != output.cols()) {
			return false;
		}

		output = *jacobian;
		return output.array().isFinite().all();
	}

	std::optional<Eigen::VectorXd> physicalStepPerSigmaFromChart(const LovCoordinateChart& chart,
	                                                             const Eigen::VectorXd& coordinates,
	                                                             const Eigen::VectorXd& coordinateStepPerSigma)
	{
		if (!chart.valid || coordinates.size() != chart.nominalCoordinates.size() ||
		    coordinateStepPerSigma.size() != coordinates.size())
			return std::nullopt;

		const auto jac = physicalJacobianFromChart(chart, coordinates);
		if (!jac)
			return std::nullopt;

		return (*jac) * coordinateStepPerSigma;
	}

}  // namespace orbital_core::hazard
