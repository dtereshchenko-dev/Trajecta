#include "orbital_core/ObservationIdentity.h"
#include <atomic>
#include <cstdint>
#include <iomanip>
#include <random>
#include <sstream>

namespace orbital_core {
	ObservationId newObservationId()
	{
		static const std::string prefix = [] {
			std::random_device random;
			std::ostringstream out;
			out << std::hex << std::setfill('0');
			for (int i = 0; i < 4; ++i)
				out << std::setw(8) << static_cast<std::uint32_t>(random());
			return out.str();
		}();
		static std::atomic<std::uint64_t> next{1};
		return prefix + '-' + std::to_string(next.fetch_add(1, std::memory_order_relaxed));
	}
}  // namespace orbital_core
