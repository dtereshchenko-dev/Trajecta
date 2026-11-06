#include "orbital_core/nongrav/NonGravProfile.h"
#include <algorithm>

NonGravProfile::NonGravProfile(NonGravModel& model) noexcept : p_model(&model) {}

NonGravModel& NonGravProfile::model() noexcept { return *p_model; }

const NonGravModel& NonGravProfile::model() const noexcept { return *p_model; }

void NonGravProfile::extendParams(std::vector<double>& params, std::vector<double>& deltaVals) const
{
	const int n = p_model->numFitParams();
	if (n == 0)
		return;

	for (int k = 0; k < n; ++k) {
		params.push_back(p_model->fitParameterValue(k));
		deltaVals.push_back(p_model->fitParameterAt(k) == NgParameter::DT ? 1.0 : 1e-14);
	}
}

void NonGravProfile::updateFromParams(NonGravModel& model, const std::vector<double>& params, int offset)
{
	const int available = static_cast<int>(params.size()) - offset;
	const int n = std::min(model.numFitParams(), std::max(0, available));
	if (n == 0)
		return;

	for (int k = 0; k < n; ++k)
		model.setFitParameterValue(k, params[offset + k]);

	if (model.model() == NgModel::Yarkovsky)
		model.setParams(0.0, model.A2(), 0.0);
}

void NonGravProfile::updateFromParams(const std::vector<double>& params, int offset) const
{
	updateFromParams(*p_model, params, offset);
}
