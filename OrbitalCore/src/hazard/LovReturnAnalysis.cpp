#include "hazard/LovReturnAnalysis.h"
#include <set>
#include <span>
#include <utility>
#include <vector>

#include <algorithm>
#include <unordered_set>

namespace orbital_core::hazard {
	namespace {
		[[nodiscard]] bool validConfig(const LovReturnAnalysisConfig& cfg) noexcept
		{
			return cfg.showerSpanDays > 0.0 && cfg.showerGapDays > 0.0;
		}

		[[nodiscard]] std::vector<int> sortByTime(std::span<const LovTraceEntry> entries)
		{
			std::vector<int> order(entries.size());
			for (int i = 0; i < static_cast<int>(entries.size()); ++i)
				order[static_cast<std::size_t>(i)] = i;

			std::sort(order.begin(), order.end(), [&](int a, int b) {
				const LovTraceEntry& ea = entries[static_cast<std::size_t>(a)];
				const LovTraceEntry& eb = entries[static_cast<std::size_t>(b)];
				if (ea.jd != eb.jd)
					return ea.jd < eb.jd;
				if (ea.sampleIndex != eb.sampleIndex)
					return ea.sampleIndex < eb.sampleIndex;
				return ea.recordIndex < eb.recordIndex;
			});
			return order;
		}

		[[nodiscard]] std::vector<std::vector<int>> splitShowers(std::span<const LovTraceEntry> entries,
		                                                         const std::vector<int>& timeOrder,
		                                                         const LovReturnAnalysisConfig& cfg)
		{
			std::vector<std::vector<int>> showers;
			if (timeOrder.empty())
				return showers;

			showers.emplace_back();
			showers.back().push_back(timeOrder.front());
			double showerStartJd = entries[static_cast<std::size_t>(timeOrder.front())].jd;
			double previousJd = showerStartJd;

			for (std::size_t i = 1; i < timeOrder.size(); ++i) {
				const int idx = timeOrder[i];
				const double jd = entries[static_cast<std::size_t>(idx)].jd;
				if (jd > showerStartJd + cfg.showerSpanDays || jd > previousJd + cfg.showerGapDays) {
					showers.emplace_back();
					showerStartJd = jd;
				}

				showers.back().push_back(idx);
				previousJd = jd;
			}
			return showers;
		}

		void sortBySampleIndex(std::span<const LovTraceEntry> entries, std::vector<int>& indices)
		{
			std::sort(indices.begin(), indices.end(), [&](int a, int b) {
				const LovTraceEntry& ea = entries[static_cast<std::size_t>(a)];
				const LovTraceEntry& eb = entries[static_cast<std::size_t>(b)];
				if (ea.sampleIndex != eb.sampleIndex)
					return ea.sampleIndex < eb.sampleIndex;
				if (ea.jd != eb.jd)
					return ea.jd < eb.jd;
				return ea.recordIndex < eb.recordIndex;
			});
		}

		void pushReturn(std::vector<LovReturnFilament>& out, const std::vector<int>& indices)
		{
			if (!indices.empty())
				out.push_back({indices});
		}

		[[nodiscard]] bool hasDuplicateSample(std::span<const LovTraceEntry> entries, const std::vector<int>& indices)
		{
			std::unordered_set<int> seen;
			seen.reserve(indices.size());
			for (const int idx : indices) {
				const int sample = entries[static_cast<std::size_t>(idx)].sampleIndex;
				if (!seen.insert(sample).second)
					return true;
			}
			return false;
		}

		void splitConsecutiveBySampleIndex(std::span<const LovTraceEntry> entries, std::vector<int> indices,
		                                   std::vector<LovReturnFilament>& out)
		{
			if (indices.empty())
				return;

			sortBySampleIndex(entries, indices);

			std::vector<int> current;
			current.push_back(indices.front());

			for (std::size_t i = 1; i < indices.size(); ++i) {
				const LovTraceEntry& prev = entries[static_cast<std::size_t>(indices[i - 1])];
				const LovTraceEntry& cur = entries[static_cast<std::size_t>(indices[i])];

				if (cur.sampleIndex > prev.sampleIndex + 1) {
					pushReturn(out, current);
					current.clear();
				}

				current.push_back(indices[i]);
			}

			pushReturn(out, current);
		}

		void splitDuplicateReturn(std::span<const LovTraceEntry> entries, std::vector<int> indices,
		                          std::vector<LovReturnFilament>& out)
		{
			if (!hasDuplicateSample(entries, indices)) {
				splitConsecutiveBySampleIndex(entries, std::move(indices), out);
				return;
			}

			std::sort(indices.begin(), indices.end(), [&](int a, int b) {
				const LovTraceEntry& ea = entries[static_cast<std::size_t>(a)];
				const LovTraceEntry& eb = entries[static_cast<std::size_t>(b)];
				if (ea.jd != eb.jd)
					return ea.jd < eb.jd;
				if (ea.sampleIndex != eb.sampleIndex)
					return ea.sampleIndex < eb.sampleIndex;
				return ea.recordIndex < eb.recordIndex;
			});

			std::vector<int> local;
			std::unordered_set<int> samples;
			samples.reserve(indices.size());

			const auto flush = [&]() {
				if (!local.empty()) {
					splitConsecutiveBySampleIndex(entries, std::move(local), out);
					local.clear();
					samples.clear();
				}
			};

			for (const int idx : indices) {
				const int sample = entries[static_cast<std::size_t>(idx)].sampleIndex;
				if (samples.contains(sample))
					flush();

				local.push_back(idx);
				samples.insert(sample);
			}
			flush();
		}

		void splitShower(std::span<const LovTraceEntry> entries, std::vector<int> indices,
		                 std::vector<LovReturnFilament>& out)
		{
			if (indices.empty())
				return;

			sortBySampleIndex(entries, indices);

			std::vector<int> current;
			current.push_back(indices.front());
			for (std::size_t i = 1; i < indices.size(); ++i) {
				const LovTraceEntry& prev = entries[static_cast<std::size_t>(indices[i - 1])];
				const LovTraceEntry& cur = entries[static_cast<std::size_t>(indices[i])];

				if (cur.sampleIndex > prev.sampleIndex + 1) {
					splitDuplicateReturn(entries, std::move(current), out);
					current.clear();
				}
				current.push_back(indices[i]);
			}

			splitDuplicateReturn(entries, std::move(current), out);
		}
	}  // namespace

	std::vector<LovReturnFilament> buildLovReturnFilaments(std::span<const LovTraceEntry> entries,
	                                                       const LovReturnAnalysisConfig& cfg)
	{
		std::vector<LovReturnFilament> out;
		if (entries.empty() || !validConfig(cfg))
			return out;

		const std::vector<int> timeOrder = sortByTime(entries);
		std::vector<std::vector<int>> showers = splitShowers(entries, timeOrder, cfg);

		for (std::vector<int>& shower : showers)
			splitShower(entries, std::move(shower), out);

		return out;
	}

}  // namespace orbital_core::hazard
