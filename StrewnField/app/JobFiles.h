#pragma once

#include "strewn_core/StrewnJob.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <string_view>

namespace strewn::jobfiles {

	using Json = nlohmann::json;

	template <std::size_t N>
	std::array<double, N> numbers(const Json& value)
	{
		if (!value.is_array() || value.size() != N)
			throw std::invalid_argument("Numeric array has the wrong dimension");
		return value.get<std::array<double, N>>();
	}

	Json readJson(const std::filesystem::path& path);
	std::string pathText(const std::filesystem::path& path);
	void writeJson(const std::filesystem::path& path, const Json& value);
	void keys(const Json& object, std::initializer_list<std::string_view> allowed);
	std::uint64_t count(const Json& object, const char* key, std::uint64_t fallback);
	EntryState readEntry(const Json& value, std::uint64_t index, bool legacy = false);
	Json entryJson(const EntryState& entry);
	const char* shapeName(Shape shape);
	Json physicsJson(const PhysicsConfig& physics);
	Json integrationJson(const IntegrationOptions& integration);
	Json samplingJson(const MonteCarloConfig& sampling);
	Json summaryJson(const RunSummary& summary);
	Environment readEnvironment(const Json& value);
	Json environmentJson(const Environment& environment);

	class ResultArchive {
	   public:
		explicit ResultArchive(const std::filesystem::path& directory);
		void append(const RunResult& result);

	   private:
		std::ofstream m_impacts;
		std::ofstream m_runs;
		std::ofstream m_breakups;
	};

}  // namespace strewn::jobfiles
