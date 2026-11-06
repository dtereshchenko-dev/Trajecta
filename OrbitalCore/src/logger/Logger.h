#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <sstream>
#include <functional>
#include <string>

enum class LogLevel { Debug, Info, Warning, Error };

class ORBITAL_CORE_API Logger {
   public:
	using Callback = std::function<void(LogLevel, const std::string&)>;

	static void setCallback(Callback cb);
	static void setLogFile(const std::string& path);
	static void closeLogFile();
	static void log(LogLevel level, const std::string& msg);

	static void setThreadMinLevel(LogLevel level) noexcept;
	static LogLevel threadMinLevel() noexcept;

	class Stream {
	   public:
		Stream(LogLevel level) : m_level(level) {}

		~Stream() { Logger::log(m_level, m_ss.str()); }

		template <typename T>
		Stream& operator<<(const T& val)
		{
			m_ss << val;
			return *this;
		}

	   private:
		LogLevel m_level;
		std::ostringstream m_ss;
	};
};

class ORBITAL_CORE_API ScopedLogSilence {
   public:
	explicit ScopedLogSilence(LogLevel newMin) noexcept : m_prev(Logger::threadMinLevel())
	{
		Logger::setThreadMinLevel(newMin);
	}

	~ScopedLogSilence() { Logger::setThreadMinLevel(m_prev); }

	ScopedLogSilence(const ScopedLogSilence&) = delete;
	ScopedLogSilence& operator=(const ScopedLogSilence&) = delete;

   private:
	LogLevel m_prev;
};

inline Logger::Stream Log(LogLevel level) { return Logger::Stream(level); }

#define LOG_DEBUG Logger::Stream(LogLevel::Debug)
#define LOG_INFO Logger::Stream(LogLevel::Info)
#define LOG_WARN Logger::Stream(LogLevel::Warning)
#define LOG_ERROR Logger::Stream(LogLevel::Error)
