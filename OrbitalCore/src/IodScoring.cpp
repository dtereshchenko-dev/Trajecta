#include "IodScoring.h"
#include "orbital_core/CometaryCoordinates.h"
#include <vector>
#include <cstdint>

#include "TwoBodyPropagation.h"
#include "obs/AstrometricSigma.h"
#include "orbital_core/Constants.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string_view>

namespace Calculations {
	double keplerResidualCost(const Vec3D& R1, const Vec3D& V1, double t1JD, const std::vector<ObsNode>& obs)
	{
		constexpr double kCostSentinel = 1.0e30;
		constexpr double kMinStateMagAU = 1.0e-12;
		double cost = 0.0;
		for (const auto& o : obs) {

			const double tau = Constants::GAUSS_K * (o.epochJdTdb() - t1JD);

			Vec3D Ri;
			Vec3D Vi;
			if (!detail::keplerPropUV(R1, V1, tau, Ri, Vi)) {
				return kCostSentinel;
			}

			const Vec3D los = Ri + o.obsVector;
			const double m = los.magnitude();

			if (m < kMinStateMagAU) {
				return kCostSentinel;
			}

			const Vec3D d = los * (1.0 / m) - o.vector;
			cost += d.dot(d);
		}

		return cost;
	}
}  // namespace Calculations

namespace Calculations::iod_scoring {
	namespace {

		constexpr double kNominalSunRadiusMeters = 695700000.0;
		constexpr double kSunRadiusAu = kNominalSunRadiusMeters / Constants::AU_IN_METERS;

		constexpr std::array<std::string_view, 30> kInclVsAScattergram = {
		    "      ):GJLKLMLKLMMNOfnxtmnklkon_D?2.09XC)&))& &&&5HJ; & & &",
		    "      4?HKMOOOPOOOOQSpy~ywyxrprpbIC=>,?bM0,).0  ,,@QSC&  &  ",
		    "      .@IKMNMNPOONPSUp|}vyz{mjqpeKNC9.>[H,00,).,00CUVD,),&&)",
		    "      ,BGHLNMPMONNOQSjz|uuwxlkqpcORB214Y>..&.&)&,<D[[H&    &",
		    "     &1=BJJKKMLLMMNOO]lojssxlrzzgXXM778]=, 0))&&,5B[[C.),&&&",
		    "      ,=DHJIKLLKKKKMIL_jgrxsit{zhWXM869V:.01,,),,,=VWA&&,. )",
		    "      ,:BEGHLKJJKHIKGEQddw}thksteOL:364O7,.1),), 2?UV?&)), )",
		    "     &,6ACFIHJHIHGONHDGXdtwoidpveNM91,4K20))&),, 0:ST> && ,)",
		    "      &3@DGFHHFGHHZcLCGKViifa_pveNJC0.:A..0)&.&&),8ST:,& )  ",
		    "     &)2=BDDGIIHFNfjKBJNJ]^]VVhm`GGD103=,,,).)&&&07RS<&&,&  ",
		    "      &7=CBFIIHGDPijIHX]R[dZQNbg]FB))0.5&)). ) )&,2QR;   )  ",
		    "      ,1<ABEGFEDFQikJMbfWZdWEJ_gXC1.&.,1&&) )& &&13OO:)&& & ",
		    "      &68>CEHGEAHMaeLL]aRTW^JJ[j]<3,)..0&) ))) ,)&1MP7 &    ",
		    "      )1;@>CBDD@EHTREJQPG^^XDF^l]<)&1)3&& &)&)&   8NO<&   & ",
		    "     &02::>@BCB>EHIHAEHE@`aTAAUcS:)&),,0&   && &) 6NO:     &",
		    "    && 17:8;:<=>AGDB=8<:A]]O=7FSJ2))&&) &&&       )MM1&    &",
		    "      &)56;:=:<?ADA;::=<5Z]J<2>JK.&.&  && && &    &EH)  )&) ",
		    "        615;::<<>@@?;5<2;QQL<55>?, ))&&& )   &   &&>@     & ",
		    "        ,53:=;;8<?=:7;846DBD51014& &&&&  &  &    &)::       ",
		    "      ) &0.1777:8>;448646=:400&)0.& ,    &      & &75,      ",
		    "      &),,142476:8<539<53374,  .22  &              ),       ",
		    "        .12300147761052,0131 ,)&.&&&)&&         &  0.&   )  ",
		    "       ) .2.4550,574,051.2 .,), )  && &   && &       &      ",
		    "     &   ,).112)205),33 &.))) & && &       &&       &       ",
		    "       &)0, 2.1&1.0,,.) 2)&.&& &,   &      &                ",
		    "        , &)2)).0,2) &.,,&   &  ,      &  &      &          ",
		    "        && .,&,.11,.),&& )&&)) & &)      &                  ",
		    "          && .0&,).1,)., ))) )& & &&       &       )        ",
		    "         ) &). &)0 ,.)),&&)&& &  &&          &              ",
		    "           ) ),&&,.&.&)&)) & &  &)&  &                      ",
		};

		[[nodiscard]] char tableValue(std::string_view row, int index) noexcept
		{
			return (index >= 0 && static_cast<std::size_t>(index) < row.size()) ? row[static_cast<std::size_t>(index)]
			                                                                    : ' ';
		}

		[[nodiscard]] double orbitLikelihoodAdjustment(double semimajorAxis, double inclinationRad,
		                                               double perihelionAu) noexcept
		{
			constexpr int kXBins = 60;
			constexpr int kYBins = 30;

			double xbin = semimajorAxis * 10.0;
			double ybin = inclinationRad * Constants::RAD_TO_DEG / 2.0;
			const int ixbin = static_cast<int>(xbin);
			const int iybin = static_cast<int>(ybin);

			double value = 0.0;
			if (ixbin >= 0 && iybin >= 0 && ixbin < kXBins - 1 && iybin < kYBins - 1 && perihelionAu > 1.1) {
				const auto row1 = kInclVsAScattergram[static_cast<std::size_t>(iybin)];
				const auto row2 = kInclVsAScattergram[static_cast<std::size_t>(iybin + 1)];
				const double p00 = static_cast<double>(tableValue(row1, ixbin) - ' ');
				const double p10 = static_cast<double>(tableValue(row1, ixbin + 1) - tableValue(row1, ixbin));
				const double p01 = static_cast<double>(tableValue(row2, ixbin) - tableValue(row1, ixbin));
				const double p11 = static_cast<double>(tableValue(row2, ixbin + 1) + tableValue(row1, ixbin) -
				                                       tableValue(row1, ixbin + 1) - tableValue(row2, ixbin));

				xbin -= static_cast<double>(ixbin);
				ybin -= static_cast<double>(iybin);
				value = p00 + p10 * xbin + p01 * ybin + p11 * xbin * ybin;
				if (perihelionAu < 1.2)
					value *= (perihelionAu - 1.1) * 10.0;
			}
			return value * 0.005;
		}

		[[nodiscard]] ClassicalElementsForScore elementsFromState(const Vec3D& R, const Vec3D& V) noexcept
		{
			const Vec3D eVec = eccentricityVec(R, V);
			const double ecc = eVec.magnitude();
			const Vec3D H = angMomentum(R, V);
			const double energy = specificEnergy(R, V, Constants::MU_IN_UNIT);

			double a = 0.0;
			if (std::abs(energy) > 1.0e-20)
				a = semiMajorAxis(Constants::MU_IN_UNIT, energy);

			double q = R.magnitude();
			if (orbital_core::cometary::usesCometaryElements(ecc)) {
				q = H.dot(H) / (1.0 + ecc);
			} else if (a != 0.0) {
				q = (ecc < 1.0) ? a * (1.0 - ecc) : std::abs(a) * (ecc - 1.0);
			}

			return {
			    .eccentricity = ecc,
			    .inclinationRad = inclination(H),
			    .semiMajorAxisAu = a,
			    .perihelionAu = q,
			};
		}

	}  // namespace

	double scoreFromElements(double weightedRms, std::size_t observationCount,
	                         const ClassicalElementsForScore& elements, bool heliocentric, bool interstellar) noexcept
	{
		if (!std::isfinite(weightedRms))
			return kBadScore;

		double score = weightedRms / 2.0;
		if (observationCount < 6)
			score *= static_cast<double>(static_cast<int>(observationCount) - 2) / 4.0;

		const double ecc = elements.eccentricity;
		const double incl = elements.inclinationRad;
		const double q = elements.perihelionAu;

		if (heliocentric) {
			score += ecc / 2.0;

			if (ecc > 1.01 && !interstellar)
				score += (ecc - 1.01) * 1000.0;

			if (incl > 0.5 && ecc < 0.8) {
				double highInclinationPenalty = (incl - 0.5) * 0.2;
				if (highInclinationPenalty > 0.2)
					highInclinationPenalty = 0.2;
				if (ecc > 0.6)
					highInclinationPenalty *= (0.8 - ecc) / 0.2;
				score += highInclinationPenalty;
			}

			score -= orbitLikelihoodAdjustment(elements.semiMajorAxisAu, incl, q);
		}

		if (ecc < 1.0 && q < kSunRadiusAu)
			score += 10000.0;

		return std::isfinite(score) ? score : kBadScore;
	}

	double weightedRmsKepler(const Vec3D& R, const Vec3D& V, double epochJD, const std::vector<ObsNode>& obs,
	                         orbital_core::AstrometricWeightProfile weightProfile) noexcept
	{
		const StateVectors state(R, V);
		if (!state.isValid() || obs.empty())
			return kBadScore;

		double sumSq = 0.0;
		int nResids = 0;
		for (const auto& o : obs) {
			const StateVectors atObs = propagateKeplerState(state, o.epochJdTdb() - epochJD);
			if (!atObs.isValid())
				return kBadScore;

			const Vec3D losEcl = atObs.first + o.obsVector;
			const Vec3D losEq(losEcl.x,
			                  losEcl.y * std::cos(Constants::EARTH_INCL) - losEcl.z * std::sin(Constants::EARTH_INCL),
			                  losEcl.y * std::sin(Constants::EARTH_INCL) + losEcl.z * std::cos(Constants::EARTH_INCL));
			const double rho = losEq.magnitude();
			if (rho <= 0.0 || !std::isfinite(rho))
				return kBadScore;

			double ra = std::atan2(losEq.y, losEq.x);
			if (ra < 0.0)
				ra += Constants::TWOPI;
			const double dec = std::asin(std::clamp(losEq.z / rho, -1.0, 1.0));

			const orbital_core::ObservationUncertainty uncertainty =
			    orbital_core::obs::lookupUncertainty(o, o.JD, weightProfile);
			double sigmaRa = uncertainty.sigmaRaRad;
			double sigmaDec = uncertainty.sigmaDecRad;
			if (sigmaRa <= 0.0 || !std::isfinite(sigmaRa))
				sigmaRa = Constants::ARCSEC_RAD;
			if (sigmaDec <= 0.0 || !std::isfinite(sigmaDec))
				sigmaDec = Constants::ARCSEC_RAD;

			const double rRa = std::cos(o.DEC) * Helpers::wrapPi(o.RA - ra) / sigmaRa;
			const double rDec = (o.DEC - dec) / sigmaDec;
			if (!std::isfinite(rRa) || !std::isfinite(rDec))
				return kBadScore;

			sumSq += rRa * rRa + rDec * rDec;
			nResids += 2;
		}

		return nResids > 0 ? std::sqrt(sumSq / static_cast<double>(nResids)) : kBadScore;
	}

	double scoreStateKepler(const Vec3D& R, const Vec3D& V, double epochJD, const std::vector<ObsNode>& obs,
	                        bool interstellar, orbital_core::AstrometricWeightProfile weightProfile) noexcept
	{
		const double weightedRms = weightedRmsKepler(R, V, epochJD, obs, weightProfile);
		if (!std::isfinite(weightedRms) || weightedRms >= kBadScore)
			return kBadScore;

		return scoreFromElements(weightedRms, obs.size(), elementsFromState(R, V), true, interstellar);
	}

	double weightedRmsFromTotalSigma(const std::vector<double>& residualsSigma,
	                                 const std::vector<std::uint8_t>& included) noexcept
	{
		double sumSq = 0.0;
		std::size_t nResiduals = 0;
		const std::size_t count = std::min(residualsSigma.size(), included.size());
		for (std::size_t i = 0; i < count; ++i) {
			if (!included[i])
				continue;

			const double rho = residualsSigma[i];
			if (std::isfinite(rho) && rho < 1.0e29) {
				sumSq += rho * rho;
				nResiduals += 2;
			}
		}

		return nResiduals > 0 ? std::sqrt(sumSq / static_cast<double>(nResiduals)) : kBadScore;
	}

	double scoreOrbit(const Orbit_p& orbit, double weightedRms, std::size_t observationCount,
	                  bool interstellar) noexcept
	{
		if (!orbit)
			return kBadScore;

		double q = orbit->getPerigee();
		if (!std::isfinite(q) && orbit->getSemiMajorAxis() != 0.0) {
			const double a = orbit->getSemiMajorAxis();
			const double e = orbit->getExcentricity();
			q = (e < 1.0) ? a * (1.0 - e) : std::abs(a) * (e - 1.0);
		}

		return scoreFromElements(weightedRms, observationCount,
		                         ClassicalElementsForScore{
		                             .eccentricity = orbit->getExcentricity(),
		                             .inclinationRad = orbit->getInclination(),
		                             .semiMajorAxisAu = orbit->getSemiMajorAxis(),
		                             .perihelionAu = q,
		                         },
		                         true, interstellar);
	}

}  // namespace Calculations::iod_scoring
