cmake_minimum_required(VERSION 3.21)
project(OrbitalCore LANGUAGES C CXX)

set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(_mc_web "${CMAKE_CURRENT_SOURCE_DIR}/src/export/web")
get_filename_component(_atmos_web "${CMAKE_CURRENT_SOURCE_DIR}/../StrewnField/output/web" ABSOLUTE)
file(READ "${_atmos_web}/viewer.html" _atmos_html)
string(REGEX MATCH "<style>[ \t\r\n]*:root[^<]*</style>" _atmos_style "${_atmos_html}")
string(REGEX REPLACE "^<style>|</style>$" "" _atmos_style "${_atmos_style}")
if(NOT _atmos_style)
    message(FATAL_ERROR "Atmospheric report stylesheet not found")
endif()
string(HEX "${_atmos_style}" _style_hex)
string(REGEX REPLACE "(..)" "0x\\1," MC_STYLE_BYTES "${_style_hex}")
file(READ "${_mc_web}/interval-mc.html" _mc_html)
foreach(_icon PRINT FIT)
    if(_icon STREQUAL "PRINT")
        set(_icon_name printer)
    else()
        set(_icon_name maximize)
    endif()
    file(READ "${_atmos_web}/vendor/lucide/${_icon_name}.svg" _svg)
    string(REPLACE "__${_icon}_ICON__" "${_svg}" _mc_html "${_mc_html}")
endforeach()
string(HEX "${_mc_html}" _html_hex)
string(REGEX REPLACE "(..)" "0x\\1," MC_HTML_BYTES "${_html_hex}")
foreach(_asset SCRIPT LEAFLET GRATICULE)
    if(_asset STREQUAL "SCRIPT")
        set(_path "${_mc_web}/interval-mc.js")
    elseif(_asset STREQUAL "LEAFLET")
        set(_path "${_atmos_web}/vendor/leaflet/leaflet.js")
    else()
        set(_path "${_atmos_web}/graticule.js")
    endif()
    file(READ "${_path}" _hex HEX)
    string(REGEX REPLACE "(..)" "0x\\1," MC_${_asset}_BYTES "${_hex}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_path}")
endforeach()
file(READ "${_atmos_web}/vendor/leaflet/leaflet.css" MC_LEAFLET_CSS)
file(READ "${_atmos_web}/vendor/leaflet/LICENSE" MC_LICENSE)
file(READ "${_atmos_web}/vendor/lucide/LICENSE" _icon_license)
string(APPEND MC_LICENSE "\nLucide\n${_icon_license}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${_mc_web}/interval-mc.html" "${_atmos_web}/viewer.html"
    "${_atmos_web}/vendor/leaflet/leaflet.css")
configure_file("${_mc_web}/IntervalMcAssets.h.in" generated/IntervalMcAssets.h @ONLY)
foreach(_asset HTML SCRIPT)
    if(_asset STREQUAL "HTML")
        set(_path "${_mc_web}/mov.html")
    else()
        set(_path "${_mc_web}/mov.js")
    endif()
    file(READ "${_path}" _hex HEX)
    string(REGEX REPLACE "(..)" "0x\\1," MOV_${_asset}_BYTES "${_hex}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_path}")
endforeach()
configure_file("${_mc_web}/MovAssets.h.in" generated/MovAssets.h @ONLY)

set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x86-windows" ${CMAKE_PREFIX_PATH})

# ---------------------------------------------------------------------------
# Numerics / perf compiler flags.
#
# The propagator / Jacobian / ephemeris paths are dominated by transcendental
# math in force evaluations and Eigen matrix ops (6x6 STM). Two
# groups of flags:
#
# 1) Applied to BOTH Debug and Release (arithmetic semantics):
#
#    /fp:precise + /fp:contract
#               — precise IEEE semantics PLUS fused multiply-add.  /fp:contract
#                 alone lifts the one transform that genuinely reduces error:
#                 `a*b + c` becomes a single vfmadd with one rounding instead
#                 of two.  Verified on this toolchain (cl 14.44): /fp:precise
#                 emits mulsd+addsd, /fp:precise /fp:contract emits vfmadd,
#                 identical to what /fp:fast produces.
#
#                 What we deliberately do NOT take from /fp:fast is
#                 reassociation.  Under /fp:fast the result becomes a function
#                 of the optimiser's inlining decisions: measured on 2026-09-11,
#                 adding /Ob3 + /LTCG alone — no source change — moved Apophis'
#                 semi-major axis by 7.8e-14 AU (~1.2 cm) and one 2020 VV
#                 impact probability by 1%.  Precise semantics decouple the
#                 numbers from build tuning, which matters more here than the
#                 ~6% /fp:fast is worth: this project chases percent-level
#                 discrepancies against JPL/ESA and needs a recompile to be a
#                 no-op on the physics.
#
#                 HISTORICAL NOTE.  This block previously claimed /fp:fast was
#                 REQUIRED, citing a 270-yr hyperbolic propagation that lost
#                 ~93% of its radial velocity under /fp:precise.  That claim
#                 does not survive arithmetic: over ~10⁶ steps, even fully
#                 coherent accumulation of double rounding bounds the error at
#                 ~10⁶ * 1.1e-16 ≈ 1e-10 relative — ten orders of magnitude
#                 short of losing 93%.  Whatever broke that run was a bug the
#                 flag masked, not a rounding effect; note the heap-corruption
#                 ABI bug documented further down this file, which is exactly
#                 the kind of defect that moves when flags change.  The run has
#                 NOT been repeated under the current code.  If you reinstate
#                 /fp:fast, reproduce that case first.
#
#    /arch:AVX2 — enables the FMA3 instruction set (vfmadd231sd etc.); without
#                 it /fp:contract has no FMA instruction to emit.  Requires
#                 Haswell (2013+).  Also fixes Eigen's packet width, which is
#                 an ABI input — see the EIGEN_MAX_ALIGN_BYTES note below.
#                 For older CPUs, downgrade to /arch:AVX (still 32-byte
#                 packets) or /arch:SSE2 (16-byte; the pinned Eigen macros
#                 keep the ABI stable either way).
#
# 2) Optimized configs only (Release + RelWithDebInfo — speed flags, no effect
#    on arithmetic semantics):
#
#    /GL, /LTCG — whole-program optimization: inline force and integration
#                 helpers across translation units, esp. Vec3D operations.
#                 Debug keeps per-TU object files for fast incremental builds.
#    /Gw        — per-variable COMDAT for globals; enables dead-data pruning.
#    /Ob3       — inline deeper than /O2's default /Ob2.  The hot paths are
#                 built from small inline Vec3D / force-model helpers, so
#                 inlining depth is what /O2 alone leaves on the table.  MSVC
#                 has no /O3 — /O2 IS its top level (cl rejects /O3 with
#                 D9002) — so /Ob3 with /GL+/LTCG is the ceiling here.
#                 Measured: /Ob3+/LTCG gave ~12% on both the Apophis fit and
#                 the 2020 VV LOV impact trace.  The root CMakeLists strips
#                 CMake's default /Ob2 so cl does not warn D9025 per file.
#
#    RelWithDebInfo ("Profile") gets the SAME speed flags on purpose: CMake
#    defaults it to /Ob1 and no /GL, which made profiling runs slower than the
#    Release they were meant to predict.  /INCREMENTAL:NO is required there
#    because CMake's RelWithDebInfo linker flags request /INCREMENTAL, which
#    /LTCG cannot honour (LNK4075).
# ---------------------------------------------------------------------------
if(MSVC)
    # Eigen/AutoDiff instantiations exceed COFF's section limit in Debug.
    set_property(SOURCE src/hazard/mov/MovCoordinates.cpp APPEND PROPERTY COMPILE_OPTIONS /bigobj)

    # C is enabled here, after the root's C++ flag cleanup.
    foreach(cfg RELEASE RELWITHDEBINFO)
        string(REGEX REPLACE "/Ob[0-9]" "" CMAKE_C_FLAGS_${cfg}
               "${CMAKE_C_FLAGS_${cfg}}")
        string(STRIP "${CMAKE_C_FLAGS_${cfg}}" CMAKE_C_FLAGS_${cfg})
    endforeach()

    # Unconditional: arithmetic semantics, both Debug and Release.
    add_compile_options(
        /fp:precise
        /fp:contract
        /arch:AVX2
    )
    # Optimized configs only: pure speed flags.
    add_compile_options(
        $<$<CONFIG:Release,RelWithDebInfo>:/GL>
        $<$<CONFIG:Release,RelWithDebInfo>:/Gw>
        $<$<CONFIG:Release,RelWithDebInfo>:/Ob3>
    )
    add_link_options(
        $<$<CONFIG:Release,RelWithDebInfo>:/LTCG>
        $<$<CONFIG:RelWithDebInfo>:/INCREMENTAL:NO>
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Same arithmetic contract as the MSVC branch above, spelled for GCC/Clang.
    #
    #  -ffp-contract=fast  is /fp:contract: FMA allowed, reassociation still
    #       forbidden (that would need -ffast-math, which this project refuses
    #       for the reasons given above).  It is NOT redundant: GCC defaults to
    #       'fast' only in gnu++ mode, and drops to 'off' under a strict
    #       -std=c++20, which CMAKE_CXX_EXTENSIONS=OFF would select -- the DLL
    #       would then quietly lose every fused multiply-add.
    #  -mavx2 -mfma  is /arch:AVX2: the FMA3 instructions the contract needs,
    #       and the same 32-byte packet width Eigen sees on Windows.  The
    #       EIGEN_MAX_*_ALIGN_BYTES pinning further down keeps the ABI fixed
    #       regardless, so a build for an older CPU only needs these two
    #       dropped.
    #
    # Deliberately absent: -flto (the /GL + /LTCG counterpart).  It changes
    # speed, not arithmetic, and costs link time on the Qt UI; enable it with
    # -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON when that trade is wanted.
    add_compile_options(
        -ffp-contract=fast
        -mavx2
        -mfma
    )
    # /Gw counterpart: per-symbol sections so the linker can drop dead data.
    add_compile_options(
        $<$<CONFIG:Release,RelWithDebInfo>:-ffunction-sections>
        $<$<CONFIG:Release,RelWithDebInfo>:-fdata-sections>
    )
    add_link_options(
        $<$<CONFIG:Release,RelWithDebInfo>:LINKER:--gc-sections>
    )
endif()

# Eigen 3.x only, and deliberately the same 3.4.x on every platform. Eigen 5
# changed reduction order and alignment; building Windows against 3.4 and Linux
# against 5.0 would reintroduce exactly the cross-platform numerical drift this
# project spends its time eliminating. Distributions that ship 5.x (Arch, for
# one) therefore need a 3.4 tree pointed at by CMAKE_PREFIX_PATH -- Eigen is
# header-only, so unpacking one is the whole installation. The linux-* presets
# already point at $ENV{HOME}/eigen-3.4.1.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR
        "Eigen 3.3+ (3.x) not found.
"
        "  A 5.x system package does NOT satisfy this: its config rejects a 3.3
"
        "  request, and mixing major versions across platforms changes results.
"
        "  Point CMAKE_PREFIX_PATH at a 3.4 tree, e.g.
"
        "      cmake --preset linux-Release -DCMAKE_PREFIX_PATH=/path/to/eigen-3.4.1")
endif()
find_package(OpenMP REQUIRED)
find_path(CSPICE_INCLUDE_DIR SpiceUsr.h PATH_SUFFIXES cspice REQUIRED)
find_library(CSPICE_LIBRARY NAMES cspice REQUIRED)
if(WIN32 AND CSPICE_LIBRARY MATCHES "\\.lib$")
    add_library(CSPICE::CSPICE SHARED IMPORTED)
    find_file(CSPICE_DLL NAMES cspice.dll HINTS "${CMAKE_PREFIX_PATH}/bin" "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin" REQUIRED)
    set_target_properties(CSPICE::CSPICE PROPERTIES IMPORTED_IMPLIB "${CSPICE_LIBRARY}" IMPORTED_LOCATION "${CSPICE_DLL}")
else()
    add_library(CSPICE::CSPICE UNKNOWN IMPORTED)
    set_target_properties(CSPICE::CSPICE PROPERTIES IMPORTED_LOCATION "${CSPICE_LIBRARY}")
endif()

# Pinned, unmodified ERFA subset for the modern IAU2006/2000A dynamical pole.
set(_erfa_root "${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa")
file(STRINGS "${_erfa_root}/sources.txt" _erfa_sources)
list(TRANSFORM _erfa_sources PREPEND "${_erfa_root}/")
add_library(OrbitalErfa STATIC ${_erfa_sources})
set_target_properties(OrbitalErfa PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES)
target_include_directories(OrbitalErfa PRIVATE "${_erfa_root}/src")

# ---------------------------------------------------------------------------
option(ORBITAL_CORE_BUILD_HAZARD
    "Build NEA hazard analysis suite (MOID, LOV, impact risk...)" ON)

# ---------------------------------------------------------------------------
# Shared library: OrbitalCore
# ---------------------------------------------------------------------------
add_library(OrbitalCore SHARED
    # ----- public API headers (3rd-party visible via #include <orbital_core/...>) -----
    api/orbital_core/OrbitalCoreExport.h
    api/orbital_core/Vec3D.h
    api/orbital_core/Constants.h
    api/orbital_core/Matrix2D.h
    api/orbital_core/Orbit.h
    api/orbital_core/ObsNode.h
    api/orbital_core/StateVectors.h
    api/orbital_core/Residuals.h
    api/orbital_core/AstrometricWeights.h
    api/orbital_core/ForceModelPolicy.h
    api/orbital_core/ForceModelSettings.h
    api/orbital_core/Types.h
    api/orbital_core/IEphemerisSource.h
    api/orbital_core/PlanetaryBodyMask.h
    api/orbital_core/DE440Source.h
    api/orbital_core/SpkReader.h

    api/orbital_core/FitSession.h
    api/orbital_core/FitRunControl.h
    api/orbital_core/ObservationSet.h
    api/orbital_core/FitPipeline.h
    api/orbital_core/FitPipelineOptions.h
    api/orbital_core/RuntimeData.h
    api/orbital_core/FitConfig.h
    api/orbital_core/FitResult.h
    api/orbital_core/Exporter.h
    api/orbital_core/forces/IForceModel.h
    api/orbital_core/forces/ForceField.h
    api/orbital_core/forces/PointMassHelpers.h
    api/orbital_core/forces/SunGravity.h
    api/orbital_core/forces/EihRelativity.h
    api/orbital_core/forces/SunJ2.h
    api/orbital_core/forces/EarthJ2.h
    api/orbital_core/perturbations/PlanetaryPerturbations.h
    api/orbital_core/perturbations/AsteroidPerturbations.h
    api/orbital_core/nongrav/NonGravModel.h
    api/orbital_core/nongrav/NonGravProfile.h
    api/orbital_core/nongrav/NonGravForceWrapper.h
    api/orbital_core/uncertainty/DerivedSigmas.h
    api/orbital_core/optimizers/ILSQProblem.h
    api/orbital_core/optimizers/EigenLevenbergMarquardt.h

    # ----- internal impl + private headers -----
    src/Types.cpp
    api/orbital_core/TimeScale.h    src/TimeScale.cpp
    src/Vec3D.cpp
    src/Orbit.cpp
    api/orbital_core/CometaryCoordinates.h
    src/CometaryCoordinates.cpp
    src/OrbitOptimizer.h            src/OrbitOptimizer.cpp
    src/Residuals.cpp
    src/EphemeridesContainer.h      src/EphemeridesContainer.cpp
    api/orbital_core/Ephemeris.h    src/ephemerides/Ephemeris.cpp
    src/obs/IObservationReader.h    src/obs/ObservationReaders.cpp
    src/obs/RwoReader.h             src/obs/RwoReader.cpp
    src/obs/Mpc80ColumnReader.h     src/obs/Mpc80ColumnReader.cpp
    api/orbital_core/obs/RadarObservation.h
    src/obs/RadarPrediction.h      src/obs/RadarPrediction.cpp
    src/AdesParser.h                src/AdesParser.cpp
    src/MPECParser.h                src/MPECParser.cpp
    src/obs/ObservationSet.cpp
    api/orbital_core/ObservationDiagnostics.h src/obs/ObservationDiagnostics.cpp
    src/obs/ObservationIdentity.cpp
    api/orbital_core/obs/AstrometricConvention.h src/obs/GaiaPhotocenter.h
    api/orbital_core/obs/GaiaFpr.h src/obs/GaiaFpr.cpp
    api/orbital_core/obs/ObservatoryUpdates.h src/obs/ObservatoryUpdates.cpp
    api/orbital_core/obs/SbdbPhysical.h src/obs/SbdbPhysical.cpp
    api/orbital_core/obs/ParsedObservation.h
    api/orbital_core/obs/StarCatalog.h src/obs/StarCatalog.cpp
    src/MPECInitData.h
    src/ObsNode.cpp
    src/FileReader.h                src/FileReader.cpp
    src/logger/Logger.h             src/logger/Logger.cpp
    src/DE440Source.cpp
    src/SpkReader.cpp
    src/io/IDataSource.h            src/io/IDataSource.cpp
    src/io/MemoryDataSource.h
    src/io/LoadedFile.h             src/io/LoadedFile.cpp
    src/io/MmappedFile.h            src/io/MmappedFile.cpp
    src/io/DataSourceFactory.h      src/io/DataSourceFactory.cpp
    src/NBodyPropagator.h           src/NBodyPropagator.cpp
    src/NBodyPropagator_SurfaceStop.cpp
    src/propagator/NBodyDynamics.h  src/propagator/NBodyDynamics.cpp
    src/propagator/EarthJ2.cpp
    src/propagator/ForceAudit.h     src/propagator/ForceAudit.cpp
    src/EnvironmentVariable.h      src/RuntimeData.cpp
    src/earth_orientation/EarthOrientation.h src/earth_orientation/EarthOrientation.cpp
    src/earth_orientation/EopData.h src/earth_orientation/C04Data.cpp
    src/earth_orientation/DeltaTData.cpp
    src/hazard/BodyFixedFrame.cpp
    src/hazard/PlanetData.cpp
    src/hazard/SurfaceGeometry.cpp
    src/integrators/IIntegrator.h
    src/integrators/BarycentricIntegrator.h src/integrators/BarycentricIntegrator.cpp
    src/integrators/RadauCore.h
    src/integrators/RA15Integrator.h   src/integrators/RA15Integrator.cpp
    src/integrators/RA15DenseOutput.h  src/integrators/RA15DenseOutput.cpp
    src/integrators/IAS15Integrator.h  src/integrators/IAS15Integrator.cpp
    src/propagator/VariationalEqs.h
    src/nongrav/NonGravModel.cpp
    src/nongrav/NonGravProfile.cpp
    src/uncertainty/DerivedSigmas.cpp
    api/orbital_core/uncertainty/CovarianceDiagnostics.h
    src/uncertainty/CovarianceDiagnostics.cpp
    src/optimizers/EigenLevenbergMarquardt.cpp
    src/obs/PreparedObservation.h
    src/obs/LightTimeSeries.h
    src/obs/AstrometricResidual.h src/obs/AstrometricResidual.cpp
    src/obs/AstrometricSigma.h     src/obs/AstrometricSigma.cpp
    src/obs/Vfcc17Rules.h          src/obs/Vfcc17Rules.cpp
    src/obs/OverObservationReweight.h src/obs/OverObservationReweight.cpp
    src/obs/ObservationPipeline.h  src/obs/ObservationPipeline.cpp
    src/obs/ObservatoryCatalog.h   src/obs/ObservatoryCatalog.cpp
    src/obs/ObserverGeometry.h     src/obs/ObserverGeometry.cpp

    src/AstrometryBias.h            src/AstrometryBias.cpp
    src/earth_orientation/EarthRotation.h src/earth_orientation/EarthRotation.cpp
    src/earth_orientation/EarthPck.h src/earth_orientation/EarthPck.cpp
    src/earth_orientation/EarthPckPoleCache.h src/earth_orientation/EarthPckPoleCache.cpp
    src/PhysicalParams.h            src/PhysicalParams.cpp
    api/orbital_core/Photometry.h    src/Photometry.cpp
    src/Gauss.h                     src/Gauss.cpp
    src/Vaisala.h                   src/Vaisala.cpp
    src/Herget.h                    src/Herget.cpp
    src/TwoBodyPropagation.h        src/TwoBodyPropagation.cpp
    src/TwoBodyStumpff.h
    src/IodScoring.h                src/IodScoring.cpp
    src/InitialOrbit.h              src/InitialOrbit.cpp
    src/PerturberSelector.h         src/PerturberSelector.cpp
    src/IODExtender.h               src/IODExtender.cpp
    src/jpl_eph/jpleph.cpp
    src/FitSession.cpp
    src/fit/SolutionRefinement.cpp
    src/fit/ObservationSelection.cpp
    src/fit/FitPassRunner.h       src/fit/FitPassRunner.cpp
    src/fit/FitPipeline.cpp        src/fit/FitPipelineDetail.h
    src/export/Exporter.cpp
    src/export/OrbitSolution.cpp
    src/fit/FitEpochPolicy.h       src/fit/FitEpochPolicy.cpp
    src/fit/HystereticRejector.h    src/fit/HystereticRejector.cpp
    src/fit/OppositionCoverage.h
    src/fit/PostFitResidual.h      src/fit/PostFitResidual.cpp
    src/fit/OptimizerConfigFactory.h src/fit/OptimizerConfigFactory.cpp
    src/fit/OrbitFitConfig.h       src/fit/OrbitFitLinearization.h
    src/fit/FitCancellation.h
    src/fit/FitObservations.h      src/fit/FitObservations.cpp
    src/fit/OrbitFitModel.h        src/fit/OrbitFitModel.cpp
    src/fit/OrbitFitProblem.h      src/fit/OrbitFitProblem.cpp
    src/fit/OrbitFitSolver.h       src/fit/OrbitFitSolver.cpp
    src/fit/OrbitFitUncertainty.h  src/fit/OrbitFitUncertainty.cpp
    src/fit/OrbitFitOutput.h       src/fit/OrbitFitOutput.cpp
    src/fit/APosterioriCovarianceScaling.h src/fit/APosterioriCovarianceScaling.cpp
)

target_sources(OrbitalCore PRIVATE
    api/orbital_core/hazard/Planet.h
    api/orbital_core/hazard/BoundaryCrossing.h
    api/orbital_core/hazard/PlanetData.h
    api/orbital_core/hazard/Moid.h
    api/orbital_core/hazard/CloseApproach.h
    api/orbital_core/hazard/ImpactEnergy.h
    api/orbital_core/hazard/PhaClassifier.h
    api/orbital_core/hazard/ImpactRisk.h
    api/orbital_core/hazard/ImpactRiskRmsObjective.h
    api/orbital_core/hazard/BPlane.h
    api/orbital_core/hazard/BodyFixedFrame.h
    api/orbital_core/hazard/SurfaceGeometry.h
    api/orbital_core/hazard/IntervalMonteCarlo.h
    api/orbital_core/hazard/LineOfVariationImpact.h
)

find_package(nlohmann_json CONFIG REQUIRED)
if(NOT TARGET TrajectaHttp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../Common/Http TrajectaHttp_build)
endif()
target_link_libraries(OrbitalCore PRIVATE nlohmann_json::nlohmann_json TrajectaHttp)

if(ORBITAL_CORE_BUILD_HAZARD)
    target_sources(OrbitalCore PRIVATE
        src/export/IntervalMcArchive.cpp
        src/export/IntervalMcReport.cpp
        src/hazard/Moid.cpp
        src/hazard/CloseApproach.cpp
        api/orbital_core/hazard/CloseApproachSearch.h
        src/hazard/CloseApproachSearch.cpp
        src/hazard/TargetPlaneDistanceBounds.h
        src/hazard/CloseApproachKinematics.h
        src/hazard/CloseApproachKinematics.cpp
        src/hazard/PhaClassifier.cpp
        src/hazard/ImpactRiskRmsObjective.cpp
        src/hazard/EncounterPlane.h
        src/hazard/EncounterPlane.cpp
        src/hazard/OpikTargetPlane.h
        src/hazard/OpikTargetPlane.cpp
        src/hazard/LineOfVariation.h
        src/hazard/LineOfVariation.cpp
        src/hazard/LovCoordinateChart.h
        src/hazard/LovCoordinateChart.cpp
        api/orbital_core/hazard/EquinoctialCoordinates.h
        src/hazard/EquinoctialCoordinates.cpp
        src/hazard/LovReturnAnalysis.h
        src/hazard/LovReturnAnalysis.cpp
        src/hazard/CloseApproachTargetPlaneDifferential.h
        src/hazard/CloseApproachTargetPlaneDifferential.cpp
        src/hazard/TargetPlaneTrace.h
        src/hazard/TargetPlaneTrace.cpp
        src/hazard/LovImpactProbability.h
        src/hazard/LovImpactProbability.cpp
        src/hazard/LineOfVariationImpact.cpp
        src/hazard/ObservationNormalMetric.h
        src/hazard/ObservationNormalMetric.cpp
        src/hazard/lov/LovTraceDetail.h
        src/hazard/lov/LovTraceGeometry.cpp
        src/hazard/lov/LovEncounterSelection.cpp
        src/hazard/lov/LovParameters.cpp
        src/hazard/lov/LovCatalog.cpp
        src/hazard/lov/LovSampling.cpp
        src/hazard/lov/LovRefinement.cpp
        src/hazard/lov/LovMinimumSearch.cpp
        src/hazard/lov/LovReturnBuilding.cpp
        src/hazard/lov/LovReturnMinima.cpp
        src/hazard/lov/LovRowAssembly.cpp
        src/hazard/VirtualImpactorRefinement.h
        src/hazard/VirtualImpactorRefinement.cpp
        src/hazard/BPlane.cpp
        src/hazard/IntervalMonteCarlo.cpp
        src/hazard/ManifoldOfVariations.cpp
    src/export/MovReport.cpp
        src/hazard/mov/MovCoordinates.cpp
        src/hazard/McBoundaryScan.cpp
        src/hazard/BoundaryScan.cpp
        src/hazard/BoundaryScan.h
        src/hazard/McObservationSampler.cpp
        src/hazard/ImpactEnergy.cpp
    )
    target_compile_definitions(OrbitalCore PUBLIC ORBITAL_CORE_HAZARD_AVAILABLE)
endif()

target_include_directories(OrbitalCore
    PRIVATE "${_erfa_root}/src" "${CMAKE_CURRENT_BINARY_DIR}/generated"
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/api
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ephemerides
        ${CMAKE_CURRENT_SOURCE_DIR}/src/jpl_eph
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io
)

target_include_directories(OrbitalCore PRIVATE "${CSPICE_INCLUDE_DIR}")
target_link_libraries(OrbitalCore PRIVATE OrbitalErfa CSPICE::CSPICE PUBLIC Eigen3::Eigen)

# Pin Eigen's ABI for every consumer of this library.
#
# /arch:AVX2 is applied by add_compile_options() above, which is DIRECTORY
# scoped: sibling targets (TrajectaCLI, TrajectaUI) do not inherit
# it.  Eigen derives EIGEN_MAX_ALIGN_BYTES from the widest available packet, so
# the DLL was built with 32 and the executables with 16.  Any type holding a
# fixed-size Eigen member whose byte size is a multiple of 16 then differs
# across the boundary -- Eigen::Matrix<double,6,6> made
# NBodyPropagator::PropagationResult 352 bytes/align 32 inside the DLL and
# 336 bytes/align 16 outside it.  A std::vector of those allocated in the DLL
# and destroyed in the caller corrupted the heap (0xC0000374), and reading
# .Phi from the caller would have silently used the wrong offset.
#
# Both macros are needed and were measured: EIGEN_MAX_ALIGN_BYTES alone leaves
# alignof(Matrix<double,6,6>) at 16 without AVX, because the alignment of
# fixed-size types comes from EIGEN_MAX_STATIC_ALIGN_BYTES.  With both pinned,
# /arch:AVX2 and plain x64 agree on 352 bytes / align 32.  Alignment wider than
# the available packet size is legal and merely costs a little padding.
target_compile_definitions(OrbitalCore PUBLIC
    EIGEN_MAX_ALIGN_BYTES=32
    EIGEN_MAX_STATIC_ALIGN_BYTES=32)

# Propagate the arithmetic flags to every consumer, which the DIRECTORY-scoped
# add_compile_options() above cannot reach.
#
# The macros pinned just above fix the ABI half of the directory-scope split.
# They do not fix the arithmetic half: much of this library is header-only --
# Vec3D's inline operators, Matrix2D, the ForceModel templates under
# api/orbital_core/forces/, the orbit factory in Orbit.h -- and those bodies
# are compiled INSIDE each consumer.  Without these flags the CLI and the UI
# evaluated them with no FMA and a different packet width than the DLL did.
#
# /fp:precise needs no propagation (it is MSVC's default) but /fp:contract
# does, otherwise consumers lose the fused multiply-add the DLL uses.
# /arch:AVX2 propagates for both FMA and Eigen's packet width.
# /GL and /LTCG are deliberately NOT propagated: whole-program optimization
# across the Qt UI buys little and costs a lot of link time.
if(MSVC)
    target_compile_options(OrbitalCore INTERFACE
        /fp:precise
        /fp:contract
        /arch:AVX2
        $<$<CONFIG:Release,RelWithDebInfo>:/Ob3>
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Same propagation argument as the MSVC branch: header-inlined bodies are
    # compiled inside each consumer, so they need the contract and the packet
    # width the library itself was built with.
    target_compile_options(OrbitalCore INTERFACE
        -ffp-contract=fast
        -mavx2
        -mfma
    )
endif()


if(OpenMP_CXX_FOUND)
    target_link_libraries(OrbitalCore PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(OrbitalCore PRIVATE ORBITAL_CORE_EXPORTS)

# ---------------------------------------------------------------------------
# Runtime data deployment
#
# orbital_core_deploy_data(<target>)
#   POST_BUILD copies DE440 kernels, MPC observatory list and the
#   reference FCCT14 and Farnocchia/Eggl 2018 debias tables next to the
#   executable.  The selected weighting profile determines which one is used.
#   Uses copy_if_different so re-builds don't touch unchanged files.
#
# Data lives in the umbrella project's data/ directory (sibling of OrbitalCore).
# Override with -DORBITAL_CORE_DATA_DIR=... if the layout differs.
# ---------------------------------------------------------------------------
set(ORBITAL_CORE_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data"
    CACHE PATH "Path to OrbitalCore runtime data directory")

function(orbital_core_deploy_data target)
    set(_d "${ORBITAL_CORE_DATA_DIR}")
    set(_dest "$<TARGET_FILE_DIR:${target}>")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_d}/linux_p1550p2650.440")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_d}/OBSCODE.dat")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_d}/RADCODE.dat")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${_d}/weights/vfcc17.rules")
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${_dest}/debias"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${_dest}/weights"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${_dest}/earth_orientation"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${_dest}/licenses"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/third_party/erfa/LICENSE" "${_dest}/licenses/ERFA.txt"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/linux_p1550p2650.440"   "${_dest}/linux_p1550p2650.440"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/OBSCODE.dat"             "${_dest}/OBSCODE.dat"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/RADCODE.dat"             "${_dest}/RADCODE.dat"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/debias/efcc18.dat"    "${_dest}/debias/efcc18.dat"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/debias/fcct14.dat"    "${_dest}/debias/fcct14.dat"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${_d}/weights/vfcc17.rules"    "${_dest}/weights/vfcc17.rules"
        COMMAND ${CMAKE_COMMAND} -E rm -f
            "${_dest}/weights/veres_2017.dat"
    )
    if(target STREQUAL "OrbitalCoreTests")
        foreach(_c04 eopc04.62-now eopc04.cnc ET-UT.dat)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${_d}/earth_orientation/${_c04}" "${_dest}/earth_orientation/${_c04}")
        endforeach()
        foreach(_pck earth_predict.bpc earth_historic.bpc earth_latest.bpc)
            if(EXISTS "${_d}/earth_orientation/${_pck}")
                add_custom_command(TARGET ${target} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        "${_d}/earth_orientation/${_pck}" "${_dest}/earth_orientation/${_pck}")
            endif()
        endforeach()
    endif()
endfunction()

# ---------------------------------------------------------------------------
# Tests (optional, skip when included as subdirectory)
# ---------------------------------------------------------------------------
option(ORBITAL_CORE_BUILD_TESTS "Build OrbitalCore tests" ON)

if(ORBITAL_CORE_BUILD_TESTS)

find_package(GTest CONFIG REQUIRED)
enable_testing()
include(GoogleTest)

add_executable(OrbitalCoreTests
    tests/TestMain.cpp
    tests/TypesTest.cpp
    tests/Vec3DTest.cpp
    tests/Matrix2DTest.cpp
    tests/OrbitTest.cpp
    tests/CometaryCoordinatesTest.cpp
    tests/OrbitFitCometaryUncertaintyTest.cpp
    src/fit/OrbitFitModel.cpp
    src/obs/RadarPrediction.cpp
    src/fit/OrbitFitProblem.cpp
    src/fit/OrbitFitUncertainty.cpp
    tests/MPECParserTest.cpp
    tests/RadarObservationTest.cpp
    tests/LightDeflectionTest.cpp
    tests/GaiaFprTest.cpp
    tests/SbdbPhysicalTest.cpp
    tests/LightTimeTest.cpp
    tests/AdesParserTest.cpp
    tests/ObservationDiagnosticsTest.cpp
    tests/StarCatalogTest.cpp
    tests/CovarianceDiagnosticsTest.cpp
    tests/PhotometryTest.cpp
    tests/Vfcc17RulesTest.cpp
    tests/OrbitOptimizerTest.cpp
    tests/GaussTest.cpp
    tests/HergetTest.cpp
    tests/VaisalaTest.cpp
    tests/EphemeridesContainerTest.cpp
    tests/EphemerisTest.cpp
    tests/DE440SourceTest.cpp
    tests/BodyFixedFrameTest.cpp
    tests/SpkReaderTest.cpp
    tests/NBodyPropagatorTest.cpp
    tests/MonteCarloAltitudeTest.cpp
    tests/IntervalMonteCarloTest.cpp
    tests/ManifoldOfVariationsTest.cpp
    src/hazard/mov/MovCoordinates.cpp
    src/hazard/McBoundaryScan.cpp
    src/hazard/BoundaryScan.cpp
    tests/IAS15IntegratorTest.cpp
    tests/RA15EncounterEpochTest.cpp
    src/integrators/RA15Integrator.cpp
    src/integrators/RA15DenseOutput.cpp
    src/integrators/IAS15Integrator.cpp
    tests/NonGravModelTest.cpp
    tests/NonGravProfileTest.cpp
    tests/EigenLevenbergMarquardtTest.cpp
    tests/AstrometryBiasTest.cpp
    tests/NutationTest.cpp
    tests/ResidualsTest.cpp
    tests/EarthJ2Test.cpp

    tests/EarthOrientationTest.cpp
    tests/HistoricalTimeTest.cpp
    tests/IDataSourceTest.cpp
    tests/FitSessionTest.cpp
    tests/ForceModelTest.cpp

    tests/OverObservationReweightTest.cpp
    tests/FitObservationsTest.cpp
    tests/APosterioriCovarianceScalingTest.cpp
    tests/EnvironmentVariableTest.cpp
    tests/ExporterOefTest.cpp
    tests/OrbitSolutionTest.cpp
    tests/ExporterRwoTest.cpp
    tests/ObservatoryUpdatesTest.cpp
    tests/PostFitResidualTest.cpp

)

if(ORBITAL_CORE_BUILD_HAZARD)
    target_sources(OrbitalCoreTests PRIVATE
        src/hazard/LineOfVariation.cpp
        src/hazard/TargetPlaneTrace.cpp
        src/hazard/ObservationNormalMetric.cpp
        src/hazard/VirtualImpactorRefinement.cpp
        src/hazard/LovReturnAnalysis.cpp
        src/hazard/lov/LovTraceGeometry.cpp
        src/hazard/lov/LovEncounterSelection.cpp
        tests/MoidTest.cpp
        tests/ImpactRiskRmsObjectiveTest.cpp
        tests/CloseApproachTest.cpp
        tests/ImpactEnergyTest.cpp
        tests/PhaClassifierTest.cpp
        tests/LovImpactProbabilityTest.cpp
        tests/LovReturnAnalysisTest.cpp
        tests/LovTraceGeometryTest.cpp
        tests/LovCatalogTest.cpp
        tests/LovEncounterSelectionTest.cpp
        src/hazard/lov/LovMinimumSearch.cpp
        tests/LovMinimumSearchTest.cpp
        src/hazard/lov/LovReturnMinima.cpp
        tests/LovReturnMinimaTest.cpp
        tests/TargetPlaneTraceTest.cpp
        tests/VirtualImpactorRefinementTest.cpp
        tests/ObservationNormalMetricTest.cpp
        tests/LovCoordinateChartTest.cpp
        tests/EquinoctialCoordinatesTest.cpp
        tests/OpikTargetPlaneTest.cpp
    )
endif()

target_link_libraries(OrbitalCoreTests PRIVATE OrbitalCore GTest::gtest CSPICE::CSPICE)
target_compile_definitions(OrbitalCoreTests PRIVATE ORBITAL_TEST_OBSERVATORY_DIR="${ORBITAL_CORE_DATA_DIR}")
orbital_core_deploy_data(OrbitalCoreTests)
if(ORBITAL_CORE_BUILD_HAZARD)
    target_link_libraries(OrbitalCoreTests PRIVATE nlohmann_json::nlohmann_json)
endif()
target_include_directories(OrbitalCoreTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src "${CSPICE_INCLUDE_DIR}")

if(WIN32)
    add_custom_command(TARGET OrbitalCoreTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:OrbitalCore>
            $<TARGET_FILE_DIR:OrbitalCoreTests>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:OrbitalCoreTests>
            $<TARGET_FILE_DIR:OrbitalCoreTests>
        COMMAND_EXPAND_LISTS
    )
endif()

gtest_discover_tests(OrbitalCoreTests)

endif() # ORBITAL_CORE_BUILD_TESTS
