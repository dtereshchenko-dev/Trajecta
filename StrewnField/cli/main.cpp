#include "strewn_core/StrewnJob.h"
#ifdef STREWN_WITH_ORBIT
#include "strewn_orbit/OrbitJob.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#endif

#include <atomic>
#include <chrono>
#include <csignal>
#include <iostream>
#include <stop_token>
#include <thread>

namespace {
	static_assert(std::atomic<bool>::is_always_lock_free);
	std::atomic<bool> interrupted{false};

	void interrupt(int) { interrupted.store(true, std::memory_order_relaxed); }

	int artifacts(const strewn::ArtifactReport& report)
	{
		for (const auto& file : report.files)
			std::cout << "artifact=" << file << '\n';

		for (const auto& warning : report.warnings)
			std::cerr << "Map warning: " << warning << '\n';

		if (!report.error.empty())
			std::cerr << "Map export: " << report.error << '\n';

		return report.status == strewn::ArtifactStatus::Cancelled
		           ? 130
		           : (report.status == strewn::ArtifactStatus::Failed ? 4 : 0);
	}
}  // namespace

#ifdef _WIN32
int wmain(int argc, wchar_t** argv)
#else
int main(int argc, char** argv)
#endif
{
	const bool render = argc == 5 && std::filesystem::path(argv[1]) == "--render";
	if ((!render && argc != 2) || (argc == 2 && std::filesystem::path(argv[1]) == "--help")) {
		std::cout << "StrewnCLI " STREWN_VERSION "\nUsage: StrewnCLI job.json\n"
		          << "       StrewnCLI --render run-directory output-directory maps.json\n"
		          << "Fixed-state or entry-sample atmospheric simulation.\n";

#ifdef STREWN_WITH_ORBIT
		std::cout << "orbit_solution: saved work directory -> orbital MC -> atmospheric simulation.\n";
#endif
		return argc == 2 ? 0 : 2;
	}

	try {
#ifdef STREWN_WITH_ORBIT
		getExecutableDir() = std::filesystem::absolute(argv[0]).parent_path();
		DE440Source ephemeris;
		SpkReader asteroids;
		setDefaultEphemerisSource(&ephemeris);
		setDefaultAsteroidSource(&asteroids);
#endif
		std::signal(SIGINT, interrupt);
		std::stop_source cancel;
		std::jthread cancellation([&](std::stop_token stop) {
			while (!stop.stop_requested()) {
				if (interrupted.load(std::memory_order_relaxed)) {
					cancel.request_stop();
					return;
				}
				std::this_thread::sleep_for(std::chrono::milliseconds(50));
			}
		});

		if (render)
			return artifacts(
			    strewn::renderSavedRun(argv[2], argv[3], strewn::loadMapOptions(argv[4]), cancel.get_token()));

		strewn::RunControl control;
		control.cancel = cancel.get_token();
		control.progress = [](const strewn::RunSummary& progress) {
			if (progress.completed == 1 || progress.completed % 10 == 0 || progress.completed == progress.requested)
				std::cout << "completed=" << progress.completed << '/' << progress.requested
				          << " impacts=" << progress.trackedImpacts << " workers=" << progress.workers
				          << " elapsed_s=" << progress.elapsedS << '\n'
				          << std::flush;

			if (!progress.convergenceHistory.empty() && progress.convergenceHistory.back().runs == progress.completed) {
				const auto& check = progress.convergenceHistory.back();
				if (check.jsCount)
					std::cout << "js_count=" << *check.jsCount << " js_mass=" << *check.jsMass
					          << " stable_checks=" << check.stableChecks << " converged=" << progress.converged << '\n'
					          << std::flush;
			}
		};

#ifdef STREWN_WITH_ORBIT
		const auto result = strewn::runConfiguredJob(std::filesystem::path(argv[1]), control, [](int done, int total) {
			if (done == 1 || done % 10 == 0 || done == total)
				std::cout << "orbital_clones=" << done << '/' << total << '\n' << std::flush;
		});
#else
		const auto result = strewn::StrewnJob{}.run(strewn::loadJob(std::filesystem::path(argv[1])), control);
#endif

		const auto mapExit = artifacts(result.artifacts);
		if (result.simulation.wallTimeouts)
			std::cerr << "wall_timeouts=" << result.simulation.wallTimeouts << '\n';

		return result.simulation.cancelled ? 130
		                                   : (result.simulation.failed || result.simulation.truncated ? 3 : mapExit);
	}

	catch (const std::exception& error) {
		std::cerr << "StrewnCLI: " << error.what() << '\n';
		return 1;
	}
}
