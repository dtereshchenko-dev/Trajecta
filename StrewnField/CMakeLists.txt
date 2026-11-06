cmake_minimum_required(VERSION 3.21)
project(StrewnField VERSION 0.2.0 LANGUAGES CXX)

option(STREWN_BUILD_TESTS "Build atmospheric flight tests" ON)
find_package(nlohmann_json CONFIG REQUIRED)
# Keep the package's debugger visualizer out of Strewn's IDE source tree.
get_target_property(_strewn_json_sources nlohmann_json::nlohmann_json INTERFACE_SOURCES)
if(_strewn_json_sources)
    list(FILTER _strewn_json_sources EXCLUDE REGEX "(^|[/\\\\])nlohmann_json[.]natvis$")
    set_property(TARGET nlohmann_json::nlohmann_json PROPERTY INTERFACE_SOURCES "${_strewn_json_sources}")
endif()
unset(_strewn_json_sources)
find_package(Threads REQUIRED)
find_package(CURL 7.85 REQUIRED)
if(NOT TARGET TrajectaHttp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../Common/Http TrajectaHttp_build)
endif()
find_package(GDAL 3.4 CONFIG REQUIRED)
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3 CONFIG REQUIRED)
endif()
include(GenerateExportHeader)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/leaflet/leaflet.js" LEAFLET_HEX HEX)
string(REGEX REPLACE "(..)" "0x\\1," LEAFLET_BYTES "${LEAFLET_HEX}")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/leaflet/leaflet.css" LEAFLET_CSS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/leaflet/LICENSE" LEAFLET_LICENSE)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/viewer.html" VIEWER_HTML)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/viewer.js" VIEWER_HEX HEX)
string(REGEX REPLACE "(..)" "0x\\1," VIEWER_BYTES "${VIEWER_HEX}")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/report.js" REPORT_JS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/print.js" PRINT_JS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/graticule.js" GRATICULE_JS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/trajectory.js" TRAJECTORY_JS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/map-flight.js" MAP_FLIGHT_JS)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/play.svg" PLAY_ICON)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/pause.svg" PAUSE_ICON)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/rotate-ccw.svg" RESTART_ICON)
foreach(ICON PLAY PAUSE RESTART)
    string(REPLACE "__${ICON}_ICON__" "${${ICON}_ICON}" VIEWER_HTML "${VIEWER_HTML}")
endforeach()
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/printer.svg" PRINT_ICON)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/maximize.svg" FIT_ICON)
string(REPLACE "__PRINT_ICON__" "${PRINT_ICON}" VIEWER_HTML "${VIEWER_HTML}")
string(REPLACE "__FIT_ICON__" "${FIT_ICON}" VIEWER_HTML "${VIEWER_HTML}")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/plotly/LICENSE" PLOTLY_LICENSE)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/lucide/LICENSE" LUCIDE_LICENSE)
string(REPLACE "__LICENSES__" "Plotly.js 3.1.0\n${PLOTLY_LICENSE}\nLucide 0.468.0\n${LUCIDE_LICENSE}" VIEWER_HTML "${VIEWER_HTML}")
string(HEX "${VIEWER_HTML}" VIEWER_HTML_HEX)
string(REGEX REPLACE "(..)" "0x\\1," VIEWER_HTML_BYTES "${VIEWER_HTML_HEX}")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/output/web/vendor/plotly/plotly.min.js" PLOTLY_HEX HEX)
string(REGEX REPLACE "(..)" "0x\\1," PLOTLY_BYTES "${PLOTLY_HEX}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    output/web/viewer.html output/web/viewer.js output/web/report.js output/web/print.js output/web/graticule.js output/web/vendor/plotly/plotly.min.js
    output/web/trajectory.js output/web/map-flight.js
    output/web/vendor/leaflet/leaflet.js output/web/vendor/leaflet/leaflet.css)
configure_file(output/web/ViewerAssets.h.in generated/ViewerAssets.h @ONLY)

add_library(StrewnCore SHARED
    api/strewn_core/AtmosphericInput.h
    api/strewn_core/Environment.h
    api/strewn_core/FieldOutput.h
    api/strewn_core/FragmentationModel.h
    api/strewn_core/Simulator.h
    api/strewn_core/StrewnJob.h
    api/strewn_core/TrajectoryIntegrator.h
    api/strewn_core/Types.h
    src/Convergence.h
    src/Random.h
    src/RunDeadline.h
    src/VectorMath.h
    src/Environment.cpp
    src/TrajectoryIntegrator.cpp
    src/FragmentationModel.cpp
    src/Simulator.cpp
    src/Convergence.cpp
    src/FieldDensity.cpp
    output/EntrySummary.h
    output/FieldOutputInternal.h
    output/ReportData.h
    output/ReportPreview.h
    output/MapOptions.cpp
    output/FieldArchive.cpp
    output/FieldOutput.cpp
    output/ImageryTiles.cpp
    output/ReportData.cpp
    output/EntrySummary.cpp
    output/MapViewer.cpp
    app/JobConfig.cpp
    app/AtmosphericInput.cpp
    app/MaterialPreset.cpp
    app/JobFiles.h
    app/JobFiles.cpp
    app/StrewnJob.cpp
    app/EnvironmentLoader.h
    app/EnvironmentLoader.cpp
    app/SourceCache.h
    app/SourceCache.cpp
    app/RasterSource.h
    app/RasterSource.cpp
    app/GfsSource.h
    app/GfsSource.cpp
    app/AutoEnvironment.h
    app/AutoEnvironment.cpp
)
target_compile_features(StrewnCore PUBLIC cxx_std_20)
set_target_properties(StrewnCore PROPERTIES CXX_EXTENSIONS OFF)
generate_export_header(StrewnCore BASE_NAME STREWN_CORE
    EXPORT_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/generated/strewn_core/Export.h")
target_include_directories(StrewnCore PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/api"
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_link_libraries(StrewnCore PRIVATE Threads::Threads nlohmann_json::nlohmann_json TrajectaHttp GDAL::GDAL
    Eigen3::Eigen)
if(WIN32)
    add_custom_command(TARGET StrewnCore POST_BUILD
        COMMAND ${CMAKE_COMMAND}
            "-DRASTER_DLL=$<TARGET_FILE:GDAL::GDAL>"
            "-DDESTINATION=$<TARGET_FILE_DIR:StrewnCore>"
            "-DCMAKE_GET_RUNTIME_DEPENDENCIES_COMMAND=${CMAKE_DUMPBIN}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeployRaster.cmake"
        VERBATIM)
endif()
if(MSVC)
    target_compile_options(StrewnCore PRIVATE /W4 /utf-8 /fp:precise)
else()
    target_compile_options(StrewnCore PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(StrewnCLI cli/main.cpp)
target_link_libraries(StrewnCLI PRIVATE StrewnCore)
target_compile_definitions(StrewnCLI PRIVATE "STREWN_VERSION=\"${PROJECT_VERSION}\"")

if(TARGET OrbitalCore AND ORBITAL_CORE_BUILD_HAZARD)
    add_library(StrewnOrbit STATIC api/strewn_orbit/OrbitJob.h orbit/OrbitJob.cpp)
    target_compile_features(StrewnOrbit PUBLIC cxx_std_20)
    target_link_libraries(StrewnOrbit PUBLIC StrewnCore OrbitalCore PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(StrewnCLI PRIVATE StrewnOrbit)
    target_compile_definitions(StrewnCLI PRIVATE STREWN_WITH_ORBIT)
    if(WIN32)
        add_custom_target(StrewnCLI-deploy-core
            COMMAND ${CMAKE_COMMAND} -E make_directory $<TARGET_FILE_DIR:StrewnCLI>
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:OrbitalCore> $<TARGET_FILE_DIR:StrewnCLI>
            DEPENDS OrbitalCore
            VERBATIM)
        add_dependencies(StrewnCLI StrewnCLI-deploy-core)
    endif()
    orbital_core_deploy_data(StrewnCLI)
endif()

if(STREWN_BUILD_TESTS)
    enable_testing()
    find_package(GTest CONFIG REQUIRED)
    include(GoogleTest)
    add_executable(StrewnCoreTests tests/CoreTest.cpp tests/JobTest.cpp tests/FieldTest.cpp tests/RandomTest.cpp tests/ConvergenceTest.cpp src/Convergence.cpp
        tests/EntrySummaryTest.cpp output/EntrySummary.cpp tests/ReportPreviewTest.cpp
        tests/HttpTest.cpp tests/WeatherTest.cpp app/GfsSource.cpp app/RasterSource.cpp app/SourceCache.cpp)
    target_include_directories(StrewnCoreTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/app")
    target_link_libraries(StrewnCoreTests PRIVATE StrewnCore GTest::gtest_main nlohmann_json::nlohmann_json TrajectaHttp GDAL::GDAL Eigen3::Eigen)
    target_compile_definitions(StrewnCoreTests PRIVATE
        STREWN_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
    gtest_discover_tests(StrewnCoreTests TEST_PREFIX Strewn. DISCOVERY_MODE PRE_TEST)
    find_program(STREWN_NODE_EXECUTABLE NAMES node)
    if(STREWN_NODE_EXECUTABLE)
        add_test(NAME Strewn.ReportStatistics
            COMMAND "${STREWN_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/ReportTest.cjs")
        add_test(NAME Strewn.ReportLoader
            COMMAND "${STREWN_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/ReportLoaderTest.cjs")
        add_test(NAME Strewn.PrintLayout
            COMMAND "${STREWN_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/PrintTest.cjs")
        add_test(NAME Strewn.Graticule
            COMMAND "${STREWN_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/GraticuleTest.cjs")
        add_test(NAME Strewn.Trajectory
            COMMAND "${STREWN_NODE_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/TrajectoryTest.cjs")
    endif()
endif()
