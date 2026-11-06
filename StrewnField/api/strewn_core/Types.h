#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <stop_token>
#include <string>
#include <variant>
#include <vector>

namespace strewn {

	using Vector3 = std::array<double, 3>;

	struct EntryState {
		Vector3 positionEcefM{};
		Vector3 velocityEcefMps{};
		double epochMjdUtc{};
		double weight{1.0};
		std::uint64_t id{};
	};

	struct FlightState {
		Vector3 positionEcefM{};
		Vector3 velocityEcefMps{};
		double massKg{};
	};

	enum class FlightPhase { Entry, DarkFlight };
	enum class Shape { Sphere, BrickAreaOnly, ConeAreaOnly };
	enum class RamPressure { RhoV2, HalfRhoV2 };
	enum class FlightEvent { Ground, Breakup, MassCutoff, Timeout, Cancelled, WallTimeout };
	enum class RunStatus { Complete, ParticleLimit, Timeout, Cancelled, Failed, WallTimeout };

	struct Fragment {
		FlightState state;
		double timeS{};
		double strengthPa{};
		double multiplicity{1.0};
		std::uint64_t id{};
		std::uint64_t parentId{};
		std::uint32_t generation{};
		Shape shape{Shape::Sphere};
	};

	struct RunParameters {
		double massKg{1.0};
		double densityKgM3{3300.0};
		double strengthPa{4e6};
		double dustFraction{0.45};
		double largestFraction{0.225};
		double beta{0.6};
		double windScale{1.0};
		double dragScale{1.0};
		double ablationKgJ{8e-8};
	};

	struct PhysicsConfig {
		FlightPhase phase{FlightPhase::Entry};
		RamPressure ramPressure{RamPressure::RhoV2};
		Shape fragmentShape{Shape::Sphere};
		double minMassKg{1e-3};
		double luminousSpeedMps{3000.0};
		double weibullAlpha{0.2};
		double strengthCapPa{330e6};
		double lateralCoefficient{1.0};
		bool dustAtEveryBreakup{true};
		std::array<std::size_t, 2> individualFragments{30, 10};
		std::array<unsigned, 2> binsPerDecade{4, 2};
		// Equal-weight directions per tail bin, capped by the bin's fragment count.
		std::size_t representativesPerBin{1};
		// Zero disables the legacy trajectory-count limit.
		std::size_t maxParticles{};
		// False reproduces the unpropagated remainder in schema-1 archives.
		bool retainMassRemainder{true};
	};

	struct IntegrationOptions {
		double relativeTolerance{1e-8};
		double positionToleranceM{0.01};
		double velocityToleranceMps{1e-4};
		double massToleranceKg{1e-9};
		double initialStepS{0.02};
		double maxStepS{5.0};
		double minStepS{1e-10};
		double maxFlightTimeS{3000.0};
		std::size_t maxStepAttempts{1000000};
	};

	struct IntegrationStats {
		std::uint64_t evaluations{};
		std::uint64_t acceptedSteps{};
		std::uint64_t rejectedSteps{};
	};

	struct FlightResult {
		FlightState state;
		double timeS{};
		FlightEvent event{FlightEvent::Ground};
		IntegrationStats stats;
	};

	struct Range {
		double min{};
		double max{};
	};

	struct MassSize {
		Range kg{1.0, 1.0};
	};

	struct DiameterSize {
		Range m{0.72, 0.72};
	};

	struct PhotometricSize {
		double absoluteMagnitude{};
		double albedoMin{0.075};
		double albedoMode{0.208};
		double albedoMax{0.33};
	};

	struct ParameterRanges {
		std::variant<MassSize, DiameterSize, PhotometricSize> size{MassSize{}};
		Range densityKgM3{2500.0, 3500.0};
		Range strongStrengthPa{0.9e6, 5e6};
		Range weakStrengthPa{0.04e6, 0.12e6};
		double weakProbability{0.0};
		Range dustFraction{0.0, 0.9};
		Range largestFraction{0.05, 0.4};
		Range beta{0.4, 0.8};
		Range windScale{0.85, 1.15};
		Range dragScale{0.85, 1.15};
		Range ablationKgJ{7e-8, 9e-8};
	};

	struct ConvergenceConfig {
		// Base-2 JS divergence of cumulative count/mass maps, not its square root.
		double threshold{0.02};
		std::size_t batchSize{50};
		std::size_t consecutiveChecks{3};
		std::array<double, 2> centerDeg{};  // Latitude, longitude; fixed for the whole run.
		std::array<double, 2> halfSizeKm{14.0, 20.0};
		double cellM{25.0};
		double kernelSigmaM{60.0};
	};

	struct MonteCarloConfig {
		std::uint64_t seed{20230213};
		std::size_t runs{200};
		unsigned workers{15};
		ParameterRanges parameters;
		std::optional<ConvergenceConfig> convergence;
		// Wall time per realization, including its entire cascade, excluding worker queue time.
		double runTimeoutS{300.0};
	};

	struct Impact {
		Fragment fragment;
		double latitudeDeg{};
		double longitudeDeg{};
		double heightWgs84M{};
	};

	struct Breakup {
		Fragment parent;
		double heightWgs84M{};
		std::size_t groups{};
	};

	struct MassBalance {
		double initialKg{};
		double landedKg{};
		double ablatedKg{};
		double dustKg{};
		double residualKg{};
		double belowCutoffKg{};
		double unresolvedKg{};
	};

	struct RunResult {
		std::uint64_t runId{};
		EntryState entry;
		RunParameters parameters;
		RunStatus status{RunStatus::Complete};
		std::size_t particles{};
		IntegrationStats integration;
		MassBalance mass;
		std::vector<Impact> impacts;
		std::vector<Breakup> breakups;
		std::string error;
	};

	struct ConvergenceCheck {
		std::size_t runs{};
		std::optional<double> jsCount;
		std::optional<double> jsMass;
		std::size_t stableChecks{};
	};

	struct RunSummary {
		std::size_t requested{};
		std::size_t completed{};
		std::size_t truncated{};
		std::size_t failed{};
		std::size_t runsWithImpacts{};
		std::size_t trackedImpacts{};
		unsigned workers{};
		bool cancelled{};
		bool converged{};
		bool convergenceBlocked{};
		std::size_t convergenceOutOfGridImpacts{};
		std::vector<ConvergenceCheck> convergenceHistory;
		double elapsedS{};
		std::size_t wallTimeouts{};
	};

	struct RunControl {
		std::stop_token cancel;
		std::function<void(const RunSummary&)> progress;
	};

	using RunConsumer = std::function<void(const RunResult&)>;
	using FlightObserver = std::function<void(double, const FlightState&)>;
	using FragmentObserver = std::function<void(const Fragment&, double, const FlightState&)>;

}  // namespace strewn
