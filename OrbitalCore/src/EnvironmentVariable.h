#pragma once

#include <cstdlib>
#include <memory>
#include <optional>
#include <string>

namespace orbital_core {
	[[nodiscard]] inline std::optional<std::string> environmentVariable(const char* name) noexcept
	{
		try {
#ifdef _WIN32
			char* raw = nullptr;
			const auto error = _dupenv_s(&raw, nullptr, name);
			const std::unique_ptr<char, decltype(&std::free)> value(raw, &std::free);
			if (error != 0 || !value)
				return std::nullopt;
			return std::string(value.get());
#else
			const char* value = std::getenv(name);
			return value ? std::optional<std::string>(value) : std::nullopt;
#endif
		} catch (...) {
			return std::nullopt;
		}
	}

#ifdef _WIN32
	[[nodiscard]] inline std::optional<std::wstring> environmentVariable(const wchar_t* name) noexcept
	{
		try {
			wchar_t* raw = nullptr;
			const auto error = _wdupenv_s(&raw, nullptr, name);
			const std::unique_ptr<wchar_t, decltype(&std::free)> value(raw, &std::free);
			if (error != 0 || !value)
				return std::nullopt;
			return std::wstring(value.get());
		} catch (...) {
			return std::nullopt;
		}
	}
#endif
}  // namespace orbital_core
