#include "EnvironmentVariable.h"

#include <gtest/gtest.h>

TEST(EnvironmentVariableTest, ReadsPresentAndAbsentValues)
{
	constexpr const char* name = "TRAJECTA_TEST_ENVIRONMENT_VALUE";
	const auto previous = orbital_core::environmentVariable(name);
	struct Restore {
		const char* name;
		std::optional<std::string> previous;
		~Restore()
		{
#ifdef _WIN32
			_putenv_s(name, previous ? previous->c_str() : "");
#else
			if (previous)
				setenv(name, previous->c_str(), 1);
			else
				unsetenv(name);
#endif
		}
	} restore{name, previous};
#ifdef _WIN32
	ASSERT_EQ(_putenv_s(name, "example"), 0);
#else
	ASSERT_EQ(setenv(name, "example", 1), 0);
#endif
	ASSERT_EQ(orbital_core::environmentVariable(name), std::optional<std::string>("example"));
#ifdef _WIN32
	ASSERT_EQ(_putenv_s(name, ""), 0);
#else
	ASSERT_EQ(unsetenv(name), 0);
#endif
	EXPECT_FALSE(orbital_core::environmentVariable(name));
}

#ifdef _WIN32
TEST(EnvironmentVariableTest, ReadsUnicodeWindowsPath)
{
	constexpr const wchar_t* name = L"TRAJECTA_TEST_ENVIRONMENT_PATH";
	const auto previous = orbital_core::environmentVariable(name);
	struct Restore {
		const wchar_t* name;
		std::optional<std::wstring> previous;
		~Restore() { _wputenv_s(name, previous ? previous->c_str() : L""); }
	} restore{name, previous};
	const std::wstring path = L"C:\\Users\\\u0422\u0435\u0441\u0442";
	ASSERT_EQ(_wputenv_s(name, path.c_str()), 0);
	EXPECT_EQ(orbital_core::environmentVariable(name), std::optional<std::wstring>(path));
}
#endif
