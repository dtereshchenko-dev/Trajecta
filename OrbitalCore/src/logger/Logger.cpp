#include "Logger.h"
#include <sstream>
#include <string>
#include <utility>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <cstdlib>

static Logger::Callback s_callback = nullptr;
static std::ofstream s_file;
static std::mutex s_mutex;

static thread_local LogLevel s_threadMinLevel = LogLevel::Info;

static const char* levelTag(LogLevel level)
{
	switch (level) {
		case LogLevel::Debug:
			return "DBG";
		case LogLevel::Info:
			return "INF";
		case LogLevel::Warning:
			return "WRN";
		case LogLevel::Error:
			return "ERR";
	}
	return "???";
}

static std::string timestamp()
{
	using namespace std::chrono;
	auto now = system_clock::now();
	auto time = system_clock::to_time_t(now);
	auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

	std::tm tm{};
#ifdef _WIN32
	localtime_s(&tm, &time);
#else
	localtime_r(&time, &tm);
#endif

	std::ostringstream ss;
	ss << std::put_time(&tm, "%H:%M:%S") << '.' << std::setfill('0') << std::setw(3) << ms.count();
	return ss.str();
}

void Logger::setCallback(Callback cb)
{
	std::lock_guard lock(s_mutex);
	s_callback = std::move(cb);
}

void Logger::setLogFile(const std::string& path)
{
	std::lock_guard lock(s_mutex);

	if (s_file.is_open())
		s_file.close();

	s_file.open(path, std::ios::trunc);

	static bool atexitRegistered = false;
	if (!atexitRegistered) {
		std::atexit(closeLogFile);
		atexitRegistered = true;
	}
}

void Logger::closeLogFile()
{
	std::lock_guard lock(s_mutex);
	if (s_file.is_open())
		s_file.close();
}

void Logger::setThreadMinLevel(LogLevel level) noexcept { s_threadMinLevel = level; }

LogLevel Logger::threadMinLevel() noexcept { return s_threadMinLevel; }

void Logger::log(LogLevel level, const std::string& msg)
{
	if (msg.empty())
		return;

	if (static_cast<int>(level) < static_cast<int>(s_threadMinLevel))
		return;

	std::lock_guard lock(s_mutex);

	if (s_callback)
		s_callback(level, msg);

	if (s_file.is_open()) {
		s_file << '[' << timestamp() << "] [" << levelTag(level) << "] " << msg << '\n';
		s_file.flush();
	}
}
