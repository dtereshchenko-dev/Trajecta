#include "strewn_core/AtmosphericInput.h"
#include <array>
#include <stdexcept>

namespace strewn {
	namespace {
		struct Material {
			Range density, strength, dust, largest, beta;
			std::array<double, 3> albedo;
		};

		constexpr std::array<Material, 3> materials{
		    {{{2500, 3500}, {.9e6, 5e6}, {0, .9}, {.05, .4}, {.4, .8}, {.075, .208, .33}},
		     {{2100, 3100}, {.25e6, .7e6}, {.3, .95}, {.01, .25}, {.5, .95}, {.02, .06, .18}},
		     {{7000, 7800}, {5e6, 50e6}, {0, .25}, {.4, .9}, {.1, .4}, {.07, .18, .54}}}};

		bool same(Range a, Range b) { return a.min == b.min && a.max == b.max; }
	}  // namespace

	void applyMaterialPreset(JobSpec& job, MaterialPreset preset)
	{
		const auto index = static_cast<std::size_t>(preset);
		if (index >= materials.size())
			throw std::invalid_argument("Choose OC, CC or Iron");
		auto* size = std::get_if<PhotometricSize>(&job.monteCarlo.parameters.size);
		if (!size)
			throw std::invalid_argument("Material preset requires the fit H for photometric size");
		const auto& material = materials[index];
		auto& p = job.monteCarlo.parameters;
		p.densityKgM3 = material.density;
		p.strongStrengthPa = material.strength;
		p.dustFraction = material.dust;
		p.largestFraction = material.largest;
		p.beta = material.beta;
		p.weakProbability = 0;
		p.weakStrengthPa = {.04e6, .12e6};
		p.windScale = p.dragScale = {.85, 1.15};
		p.ablationKgJ = {7e-8, 9e-8};
		size->albedoMin = material.albedo[0];
		size->albedoMode = material.albedo[1];
		size->albedoMax = material.albedo[2];
	}

	MaterialPreset materialPreset(const JobSpec& job)
	{
		const auto& p = job.monteCarlo.parameters;
		const auto* size = std::get_if<PhotometricSize>(&p.size);
		if (!size || p.weakProbability != 0 || !same(p.windScale, {.85, 1.15}) || !same(p.dragScale, {.85, 1.15}) ||
		    !same(p.ablationKgJ, {7e-8, 9e-8}))
			return MaterialPreset::Custom;
		for (std::size_t i = 0; i < materials.size(); ++i) {
			const auto& m = materials[i];
			if (same(p.densityKgM3, m.density) && same(p.strongStrengthPa, m.strength) &&
			    same(p.dustFraction, m.dust) && same(p.largestFraction, m.largest) && same(p.beta, m.beta) &&
			    size->albedoMin == m.albedo[0] && size->albedoMode == m.albedo[1] && size->albedoMax == m.albedo[2])
				return static_cast<MaterialPreset>(i);
		}
		return MaterialPreset::Custom;
	}
}  // namespace strewn
