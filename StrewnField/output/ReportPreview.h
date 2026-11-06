#pragma once

#include "FieldOutputInternal.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <span>

namespace strewn::output {

	// Display budgets do not enter propagation, KDE or ensemble statistics.
	inline std::size_t previewLimit(std::size_t runs)
	{
		return std::clamp(std::size_t{32000} / std::max(std::size_t{1}, runs), std::size_t{32}, std::size_t{500});
	}

	struct PreviewTrack {
		std::uint64_t id, parent;
		double mass;
	};

	// Retain complete ancestry; rank-stratified terminal masses are not a statistical sample.
	inline std::set<std::uint64_t> previewTracks(std::span<const PreviewTrack> tracks, std::size_t limit)
	{
		std::map<std::uint64_t, const PreviewTrack*> byId;
		std::set<std::uint64_t> parents, chosen;
		for (const auto& t : tracks) {
			byId.emplace(t.id, &t);
			if (t.id != t.parent)
				parents.insert(t.parent);
		}
		if (tracks.size() <= limit) {
			for (const auto& t : tracks)
				chosen.insert(t.id);
			return chosen;
		}
		for (const auto& t : tracks)
			if ((t.id == t.parent || !byId.contains(t.parent)) && chosen.size() < limit)
				chosen.insert(t.id);
		std::vector<const PreviewTrack*> leaves;
		for (const auto& t : tracks)
			if (!parents.contains(t.id))
				leaves.push_back(&t);
		std::sort(leaves.begin(), leaves.end(),
		          [](auto a, auto b) { return a->mass != b->mass ? a->mass > b->mass : a->id < b->id; });
		const auto add = [&](const PreviewTrack* node) {
			std::set<std::uint64_t> chain;
			while (node && !chosen.contains(node->id)) {
				if (!chain.insert(node->id).second)
					break;
				const auto it = byId.find(node->parent);
				node = it == byId.end() ? nullptr : it->second;
			}
			if (chosen.size() + chain.size() <= limit)
				chosen.insert(chain.begin(), chain.end());
		};
		if (!leaves.empty()) {
			add(leaves.front());
			add(leaves.back());
		}
		std::deque<std::pair<std::size_t, std::size_t>> spans;
		if (leaves.size() > 2)
			spans.emplace_back(1, leaves.size() - 2);
		while (!spans.empty() && chosen.size() < limit) {
			const auto [lo, hi] = spans.front();
			spans.pop_front();
			const auto mid = lo + (hi - lo) / 2;
			add(leaves[mid]);
			if (mid > lo)
				spans.emplace_back(lo, mid - 1);
			if (mid < hi)
				spans.emplace_back(mid + 1, hi);
		}
		return chosen;
	}

	// Identical bin definitions to StrewnReport.histogram, evaluated on the full archive.
	template <class Value>
	Json impactHistogram(const ArchivedField& archive, Value value, bool logarithmic)
	{
		constexpr std::size_t bins = 28;
		Json result = {{"x", Json::array()}, {"y", Json::array()}, {"width", Json::array()}, {"ranges", Json::array()}};
		double lo = std::numeric_limits<double>::infinity(), hi = -lo;
		for (std::size_t i = 0; i < archive.points.size(); ++i) {
			const auto& p = archive.points[i];
			const double v = value(i);
			if (p.multiplicity * p.runWeight <= 0 || !std::isfinite(v) || (logarithmic && v <= 0))
				continue;
			const double x = logarithmic ? std::log10(v) : v;
			lo = std::min(lo, x);
			hi = std::max(hi, x);
		}
		if (!std::isfinite(lo))
			return result;
		if (hi - lo <= 32 * std::numeric_limits<double>::epsilon() * std::max({1., std::abs(lo), std::abs(hi)})) {
			const double v = logarithmic ? std::pow(10, (lo + hi) / 2) : (lo + hi) / 2;
			return {{"x", {v}}, {"y", {1.}}, {"width", {0.}}, {"ranges", {{v, v}}}, {"singleValue", v}};
		}
		const double step = (hi - lo) / bins;
		std::array<double, bins> weights{};
		for (std::size_t i = 0; i < archive.points.size(); ++i) {
			const auto& p = archive.points[i];
			const double v = value(i), w = p.multiplicity * p.runWeight;
			if (w <= 0 || !std::isfinite(v) || (logarithmic && v <= 0))
				continue;
			const double x = logarithmic ? std::log10(v) : v;
			weights[std::min(bins - 1, static_cast<std::size_t>(std::max(0., std::floor((x - lo) / step))))] += w;
		}
		double total = 0;
		for (double w : weights)
			total += w;
		for (std::size_t i = 0; i < bins; ++i) {
			const double a = logarithmic ? std::pow(10, lo + i * step) : lo + i * step;
			const double b = logarithmic ? std::pow(10, lo + (i + 1) * step) : lo + (i + 1) * step;
			result["x"].push_back((a + b) / 2);
			result["width"].push_back((b - a) * .96);
			result["ranges"].push_back({a, b});
			result["y"].push_back(weights[i] / total);
		}
		return result;
	}
}  // namespace strewn::output
