#include "orbital_core/Photometry.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace orbital_core::photometry {
	namespace {

		std::string_view trim(std::string_view text) noexcept
		{
			const auto first = text.find_first_not_of(" \t\r\n");
			if (first == text.npos)
				return {};
			return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
		}

		int decimalPlaces(std::string_view text) noexcept
		{
			text = trim(text);
			const auto dot = text.find('.');
			if (dot == text.npos)
				return 0;
			int count = 0;
			for (auto i = dot + 1; i < text.size() && text[i] >= '0' && text[i] <= '9'; ++i)
				++count;
			return count;
		}

		bool positiveFinite(double value) noexcept { return std::isfinite(value) && value > 0.0; }

		std::optional<double> bandToVCorrection(char band) noexcept
		{
			struct BandShift {
				char band;
				double correction;
			};

			static constexpr std::array<BandShift, 22> shifts{{
			    {' ', -0.80}, {'B', -0.80}, {'c', -0.05}, {'C', 0.40}, {'g', -0.28}, {'G', 0.28},
			    {'H', 0.00},  {'i', 0.39},  {'I', 0.80},  {'J', 1.10}, {'K', 1.70},  {'L', 0.20},
			    {'o', 0.33},  {'r', 0.23},  {'R', 0.40},  {'u', 2.50}, {'U', 1.10},  {'w', 0.16},
			    {'y', 0.36},  {'Y', 0.70},  {'z', 0.37},  {'V', 0.00},
			}};

			for (const auto& shift : shifts) {
				if (shift.band == band)
					return shift.correction;
			}

			return std::nullopt;
		}

		std::optional<double> individualMagnitude(const ::Residual& obs, double slopeG)
		{
			if (!std::isfinite(obs.mag) || obs.mag <= 0.1 || !std::isfinite(obs.helioDistAu) ||
			    obs.helioDistAu < 1.0e-6 || !std::isfinite(obs.geoDistAu) || obs.geoDistAu < 1.0e-6 ||
			    !std::isfinite(obs.phaseAngleRad))
				return std::nullopt;

			const auto bandCorrection = bandToVCorrection(obs.magBand);
			if (!bandCorrection)
				return std::nullopt;

			double phaseCorrection = 0.0;
			if (obs.phaseAngleRad > 0.0) {
				phaseCorrection = 2.5 * std::log10(Calculations::hgPhaseFunction(obs.phaseAngleRad, slopeG));
			}

			const double h =
			    obs.mag - 5.0 * std::log10(obs.helioDistAu * obs.geoDistAu) + phaseCorrection + *bandCorrection;
			return std::isfinite(h) ? std::optional<double>{h} : std::nullopt;
		}

		double effectiveSigma(const ::Residual& obs) noexcept
		{
			return positiveFinite(obs.magSigma) ? obs.magSigma : 0.7;
		}

	}  // namespace

	double resolveSigma(std::string_view magnitudeText, double reportedSigma) noexcept
	{
		if (positiveFinite(reportedSigma))
			return reportedSigma;

		if (trim(magnitudeText).empty())
			return 0.0;

		const int digits = decimalPlaces(magnitudeText);
		return digits == 0 ? 1.0 : digits == 1 ? 0.7 : 0.5;
	}

	const Residuals& observations(const FitResult& result) noexcept
	{
		return result.allResiduals.count() != 0 ? result.allResiduals : result.residuals;
	}

	Solution fit(std::span<const ::Residual> observations, double slopeG)
	{
		Solution result;
		result.slopeG = slopeG;
		if (!std::isfinite(slopeG))
			return result;
		double hsum = 0.0;
		for (const auto& obs : observations) {
			const auto h = individualMagnitude(obs, slopeG);
			if (!h)
				continue;

			const double weight = 1.0 / effectiveSigma(obs);
			hsum += *h * weight;
			result.weightSum += weight;
			++result.count;
		}

		if (result.count == 0)
			return result;

		result.absoluteMagnitude = hsum / result.weightSum;
		double residualSum = 0.0;
		double squaredWeightSum = 0.0;
		for (const auto& obs : observations) {
			const auto h = individualMagnitude(obs, slopeG);
			if (!h)
				continue;

			const double weight = 1.0 / effectiveSigma(obs);
			const double residual = (*h - *result.absoluteMagnitude) * weight;

			residualSum += residual * residual;
			squaredWeightSum += weight * weight;
		}

		result.rms = std::sqrt(residualSum / squaredWeightSum);
		return result;
	}

	CometSolution fitComet(std::span<const ::Residual> observations, double slope)
	{
		PhotometryOptions options;
		options.cometK1 = options.cometK2 = slope;
		return fitComet(observations, options);
	}

	CometSolution fitComet(std::span<const ::Residual> observations, const PhotometryOptions& options)
	{
		CometSolution result;
		result.total.slope = options.cometK1;
		result.nuclear.slope = options.cometK2;
		result.nuclear.phaseCoefficient = options.cometNuclearPhaseCoefficient;
		if (!std::isfinite(options.cometK1) || !std::isfinite(options.cometK2) ||
		    !std::isfinite(options.cometNuclearPhaseCoefficient))
			return result;
		double sums[2]{}, weights[2]{};
		for (const auto& obs : observations) {
			if (!std::isfinite(obs.mag) || !positiveFinite(obs.magSigma) || !positiveFinite(obs.helioDistAu) ||
			    !positiveFinite(obs.geoDistAu))
				continue;
			if (!obs.nuclearMagnitude) {
				++result.unknownTypeCount;
				continue;
			}
			const auto correction = (obs.magBand == 'N' || obs.magBand == 'T') ? std::optional<double>{0.0}
			                                                                   : bandToVCorrection(obs.magBand);
			if (!correction)
				continue;
			const int index = *obs.nuclearMagnitude ? 1 : 0;
			const auto& component = index ? result.nuclear : result.total;
			const auto geometry = cometMagnitude(0.0, component.slope, obs.helioDistAu, obs.geoDistAu,
			                                     obs.phaseAngleRad, component.phaseCoefficient);
			if (!geometry)
				continue;
			const auto reduced = obs.mag + *correction - *geometry;
			const double weight = 1.0 / obs.magSigma;
			sums[index] += reduced * weight;
			weights[index] += weight;
			++(index ? result.nuclear.count : result.total.count);
		}
		if (weights[0] > 0.0)
			result.total.absoluteMagnitude = sums[0] / weights[0];
		if (weights[1] > 0.0)
			result.nuclear.absoluteMagnitude = sums[1] / weights[1];
		return result;
	}

	std::optional<double> cometMagnitude(double absoluteMagnitude, double slope, double heliocentricAu,
	                                     double observerAu, double phaseRad, double phaseCoefficient) noexcept
	{
		if (!std::isfinite(absoluteMagnitude) || !std::isfinite(slope) || !positiveFinite(heliocentricAu) ||
		    !positiveFinite(observerAu) || !std::isfinite(phaseCoefficient) ||
		    (phaseCoefficient != 0.0 && !std::isfinite(phaseRad)))
			return std::nullopt;
		const double phase = phaseCoefficient == 0.0 ? 0.0 : phaseCoefficient * phaseRad * 180.0 / std::acos(-1.0);
		return absoluteMagnitude + 5.0 * std::log10(observerAu) + slope * std::log10(heliocentricAu) + phase;
	}

	Evaluation evaluateComet(const ::Residual& obs, const CometSolution& solution, bool includedInFit)
	{
		Evaluation result;
		if (!obs.nuclearMagnitude || !std::isfinite(obs.mag) || !positiveFinite(obs.magSigma))
			return result;
		const auto& component = *obs.nuclearMagnitude ? solution.nuclear : solution.total;
		if (!component.absoluteMagnitude)
			return result;
		const auto correction =
		    (obs.magBand == 'N' || obs.magBand == 'T') ? std::optional<double>{0.0} : bandToVCorrection(obs.magBand);
		const auto predicted = cometMagnitude(*component.absoluteMagnitude, component.slope, obs.helioDistAu,
		                                      obs.geoDistAu, obs.phaseAngleRad, component.phaseCoefficient);
		if (!correction || !predicted)
			return result;
		result.sigma = obs.magSigma;
		result.residual = obs.mag + *correction - *predicted;
		result.selected = includedInFit;
		return result;
	}

	Evaluation evaluate(const ::Residual& obs, const Solution& solution, bool includedInFit)
	{
		Evaluation result;
		const auto h = individualMagnitude(obs, solution.slopeG);
		if (!h)
			return result;

		result.sigma = effectiveSigma(obs);
		if (!solution.absoluteMagnitude)
			return result;

		result.residual = *h - *solution.absoluteMagnitude;
		result.selected = true;

		if (solution.count > 4) {
			double checkResidual = *result.residual;

			if (includedInFit) {
				const double remainingWeight = solution.weightSum - 1.0 / result.sigma;
				if (remainingWeight > 0.0)
					checkResidual *= solution.weightSum / remainingWeight;
			}

			result.selected = std::pow(checkResidual / result.sigma, 2.0) <= 8.0;
		}

		return result;
	}

	std::string formatMagnitude(double magnitude, char band, std::string_view magnitudeText)
	{
		if (!positiveFinite(magnitude))
			return std::string(6, ' ');

		const int precision = trim(magnitudeText).empty() ? 2 : std::min(decimalPlaces(magnitudeText), 2);
		for (int digits = precision; digits >= 0; --digits) {
			std::ostringstream out;
			out.imbue(std::locale::classic());
			out << std::fixed << std::setprecision(digits) << std::setw(digits == 0 ? 2 : 3 + digits) << magnitude;
			std::string text = out.str();
			if (text.size() > 5)
				continue;

			text.resize(5, ' ');
			text.push_back(band == '\0' ? ' ' : band);
			return text;
		}

		return std::string(6, ' ');
	}

	std::optional<double> visualMagnitude(double absoluteMagnitude, double slopeG, double heliocentricAu,
	                                      double observerAu, double phaseRad) noexcept
	{
		if (!std::isfinite(absoluteMagnitude) || !std::isfinite(slopeG) || !positiveFinite(heliocentricAu) ||
		    !positiveFinite(observerAu) || !std::isfinite(phaseRad) || phaseRad < 0.0 ||
		    phaseRad >= 3.14159265358979323846)
			return std::nullopt;

		const double phaseFunction = Calculations::hgPhaseFunction(phaseRad, slopeG);
		if (!(phaseFunction > 0.0))
			return std::nullopt;
		return absoluteMagnitude + 5.0 * std::log10(heliocentricAu * observerAu) - 2.5 * std::log10(phaseFunction);
	}

}  // namespace orbital_core::photometry
