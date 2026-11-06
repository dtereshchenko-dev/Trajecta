#include "obs/IObservationReader.h"
#include <string_view>
#include <vector>

#include "AdesParser.h"
#include "obs/Mpc80ColumnReader.h"
#include "obs/RwoReader.h"

namespace orbital_core::obs {

	std::vector<const IObservationReader*> observationReaders() noexcept
	{
		static const ades::AdesReader ades;
		static const Mpc80ColumnReader mpc80;
		static const RwoReader rwo;
		static const std::vector<const IObservationReader*> readers = {
		    &ades,
		    &rwo,
		    &mpc80,
		};
		return readers;
	}

	const IObservationReader* selectObservationReader(std::string_view text) noexcept
	{
		for (const IObservationReader* reader : observationReaders()) {
			if (reader->canRead(text))
				return reader;
		}
		return nullptr;
	}

}  // namespace orbital_core::obs
