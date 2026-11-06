#include "EphemeridesContainer.h"
#include <memory>
#include <string>
#include "orbital_core/Types.h"
#include "orbital_core/IEphemerisSource.h"

struct EphemeridesContainer::Impl {
	Error error = Error::NoError;
	Planet planetId;

	CelestianEphemerisData getFromDE440(double target, bool isMjd = false);
};

EphemeridesContainer::EphemeridesContainer(const Planet planetId) : d(std::make_unique<Impl>())
{
	d->planetId = planetId;
}

EphemeridesContainer::~EphemeridesContainer() = default;

EphemeridesContainer::CelestianEphemerisData EphemeridesContainer::get(const double& targetJD)
{
	d->error = Error::NoError;

	auto* source = getDefaultEphemerisSource();
	if (source && source->isAvailable())
		return d->getFromDE440(targetJD);

	d->error = Error::IOError;
	return {};
}

EphemeridesContainer::Error EphemeridesContainer::error() const noexcept { return d->error; }

EphemeridesContainer::CelestianEphemerisData EphemeridesContainer::getMjd(double targetMjd)
{
	d->error = Error::NoError;
	auto* source = getDefaultEphemerisSource();
	if (source && source->isAvailable())
		return d->getFromDE440(targetMjd, true);
	d->error = Error::IOError;
	return {};
}

std::string EphemeridesContainer::errorString() const noexcept
{
	switch (d->error) {
		case Error::NoError:
			return "No error";
		case Error::InvalidData:
			return "Error: Invalid data";
		case Error::InsufficientData:
			return "Error: Insufficient data for interpolation";
		case Error::IOError:
			return "Error: Error reading file";
	}
	return "Unknown error";
}

EphemeridesContainer::CelestianEphemerisData EphemeridesContainer::Impl::getFromDE440(double target, bool isMjd)
{
	auto* source = getDefaultEphemerisSource();

	EphemerisPlanet ep;
	switch (planetId) {
		case Planet::Sun:
			ep = EphemerisPlanet::Sun;
			break;
		case Planet::Mercury:
			ep = EphemerisPlanet::Mercury;
			break;
		case Planet::Venus:
			ep = EphemerisPlanet::Venus;
			break;
		case Planet::Earth:
			ep = EphemerisPlanet::Earth;
			break;
		case Planet::Mars:
			ep = EphemerisPlanet::Mars;
			break;
		case Planet::Jupiter:
			ep = EphemerisPlanet::Jupiter;
			break;
		case Planet::Saturn:
			ep = EphemerisPlanet::Saturn;
			break;
		case Planet::Uranus:
			ep = EphemerisPlanet::Uranus;
			break;
		case Planet::Neptune:
			ep = EphemerisPlanet::Neptune;
			break;
		case Planet::Moon:
			ep = EphemerisPlanet::Moon;
			break;
		default:
			error = Error::InvalidData;
			return {};
	}

	auto rec = isMjd ? source->queryMjd(ep, target) : source->queryJD(ep, target);

	CelestianEphemerisData data{};
	data.jd = rec.jd;
	data.R = rec.R;
	data.V = rec.V;
	data.dobs = rec.dobs;
	data.phaseAngle = rec.phaseAngle;
	data.mv = rec.mv;
	data.vr = rec.vr;
	return data;
}
