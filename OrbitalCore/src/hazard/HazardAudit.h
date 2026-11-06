#pragma once

#include "logger/Logger.h"
#include "EnvironmentVariable.h"

#include <iomanip>
#include <locale>
#include <sstream>

namespace orbital_core::hazard::audit {
	struct Context {
		int sample;
		double sigma;
	};

	inline thread_local const Context* current = nullptr;

	inline bool enabled()
	{
		static const bool value = [] {
			const auto setting = environmentVariable("ORBITAL_HAZARD_AUDIT");
			return setting && !setting->empty() && (*setting)[0] == '1';
		}();
		return value;
	}

	class Scope final {
	   public:
		Scope(int sample, double sigma, bool force = false) : context{sample, sigma}, previous(current)
		{
			if (force || previous || enabled())
				current = &context;
		}

		~Scope() { current = previous; }

		Scope(const Scope&) = delete;
		Scope& operator=(const Scope&) = delete;

	   private:
		Context context;
		const Context* previous;
	};

	// One logger record per stage; sample/sigma join parallel traces.
	template <class Writer>
	void write(const char* stage, Writer&& writer)
	{
		if (!current)
			return;
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << std::setprecision(17) << "[HAZARD_AUDIT] stage=" << stage << " sample=" << current->sample
		    << " sigma=" << current->sigma;
		writer(out);
		LOG_INFO << out.str();
	}
}  // namespace orbital_core::hazard::audit
