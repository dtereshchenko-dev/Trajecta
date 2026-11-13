cmake_minimum_required(VERSION 3.21)
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/TrajectaVersion.cmake)
project(TrajectaUI VERSION ${TRAJECTA_VERSION} LANGUAGES CXX)

set(CMAKE_PREFIX_PATH "C:/Qt/6.9.3/msvc2022_64" CACHE PATH "Qt6 installation")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OrbitalCore shared library (without tests)
if(NOT TARGET OrbitalCore)
    set(ORBITAL_CORE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../OrbitalCore OrbitalCore_build)
endif()

find_package(Qt6 REQUIRED COMPONENTS
    Quick
    QuickControls2
    Quick3D
    Charts
    Widgets
    Concurrent
    Network
)

qt_standard_project_setup(REQUIRES 6.5)
if(Qt6_VERSION VERSION_GREATER_EQUAL 6.8)
    qt_policy(SET QTP0004 NEW)
endif()

list(APPEND QML_IMPORT_PATH "${CMAKE_CURRENT_BINARY_DIR}")
list(REMOVE_DUPLICATES QML_IMPORT_PATH)
set(QML_IMPORT_PATH "${QML_IMPORT_PATH}" CACHE STRING "QML module paths for Qt Creator" FORCE)

if(NOT TARGET StrewnCore)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../StrewnField StrewnField_build)
endif()

qt_add_executable(TrajectaUI
    main.cpp
    app.rc
)

if(WIN32)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/app-version.rc.in
                   ${CMAKE_CURRENT_BINARY_DIR}/app-version.rc @ONLY)
    target_sources(TrajectaUI PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/app-version.rc)
endif()
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/trajecta-version.txt CONTENT "${TRAJECTA_VERSION}\n")

set_target_properties(TrajectaUI PROPERTIES
    WIN32_EXECUTABLE TRUE
    OUTPUT_NAME Trajecta
)

set_source_files_properties(qml/shared/BodyPalette.qml PROPERTIES
    QT_QML_SINGLETON_TYPE TRUE
)

set_source_files_properties(qml/shared/AppTheme.qml qml/shared/DialogState.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)
set(CLIENT_SOURCES
    src/controllers/OrbitEngine.h src/controllers/OrbitEngine.cpp
    src/controllers/OrbitWorker.h src/controllers/OrbitWorker.cpp
    src/geometry/OrbitGeometry.h src/geometry/OrbitGeometry.cpp
    src/geometry/GridGeometry.h src/geometry/GridGeometry.cpp
    src/io/RecentFiles.h src/io/RecentFiles.cpp
    src/util/FormatUtils.h src/util/CalculationTime.h
    src/util/PlanetPositionsProvider.h src/util/PlanetPositionsProvider.cpp
    src/models/DiagnosticsModel.h src/models/DiagnosticsModel.cpp
    src/models/OrbitModel.h src/models/OrbitModel.cpp
    src/settings/PerturbationsSettings.h src/settings/PerturbationsSettings.cpp
    src/settings/NonGravSettings.h src/settings/NonGravSettings.cpp
    src/settings/WeightingSettings.h src/settings/WeightingSettings.cpp
    src/settings/PhotometrySettings.h src/settings/PhotometrySettings.cpp
    src/settings/SigmaClipSettings.h src/settings/SigmaClipSettings.cpp
    src/settings/FitSettingsAggregator.h src/settings/FitSettingsAggregator.cpp
    src/settings/SettingsStorage.h src/settings/SettingsStorage.cpp
    src/features/risk/RiskMonitorViewModel.h src/features/risk/RiskMonitorViewModel.cpp
    src/features/risk/MovViewModel.h src/features/risk/MovViewModel.cpp
    src/features/atmosphere/AtmosphereViewModel.h src/features/atmosphere/AtmosphereViewModel.cpp
    src/features/atmosphere/AtmosphereArchive.h src/features/atmosphere/AtmosphereArchive.cpp
    src/features/atmosphere/AtmosphereExecution.cpp
    src/features/risk/RiskEventTableModel.h src/features/risk/RiskEventTableModel.cpp
    src/features/solution/ResidualMarker.h src/features/solution/ResidualMarker.cpp
    src/features/orbit/OrbitSceneModel.h src/features/orbit/OrbitSceneModel.cpp
    src/features/settings/CalculationSettingsModel.h src/features/settings/CalculationSettingsModel.cpp
    src/features/settings/DataSourcesModel.h src/features/settings/DataSourcesModel.cpp
    src/features/settings/Sb441Download.h src/features/settings/Sb441Download.cpp
    src/features/settings/ObservatoryUpdate.h src/features/settings/ObservatoryUpdate.cpp
    src/features/settings/ApplicationUpdate.h src/features/settings/ApplicationUpdate.cpp
    src/app/ApplicationSession.h src/app/ApplicationSession.cpp
    src/features/solution/SolutionViewModel.h src/features/solution/SolutionViewModel.cpp
    src/features/solution/CovarianceTableModel.h src/features/solution/CovarianceTableModel.cpp
    src/features/observations/ObservationTableModel.h src/features/observations/ObservationTableModel.cpp
    src/features/observations/ObservationReviewModel.h src/features/observations/ObservationReviewModel.cpp
    src/features/observations/StationResidualModel.h src/features/observations/StationResidualModel.cpp
    src/features/ephemerides/EphemerisTableModel.h src/features/ephemerides/EphemerisTableModel.cpp
    src/features/ephemerides/EphemerisViewModel.h src/features/ephemerides/EphemerisViewModel.cpp
    src/features/encounters/CloseApproachesViewModel.h src/features/encounters/CloseApproachesViewModel.cpp
    src/features/encounters/CloseApproachTableModel.h src/features/encounters/CloseApproachTableModel.cpp
)
set(CLIENT_QML_FILES
    qml/shared/BodyPalette.qml
    qml/features/encounters/CloseApproachesPage.qml
    qml/features/encounters/CloseApproachTable.qml
    qml/features/risk/MovCalculationForm.qml
    qml/features/risk/MovRunView.qml
    qml/features/risk/MovSourcePanel.qml
    qml/components/shared/FormField.qml
    qml/components/shared/FormSection.qml
    qml/components/shared/SearchIntervalForm.qml
    qml/components/shared/EmptyResultTable.qml
    qml/features/risk/RiskMonitorPage.qml
    qml/features/risk/RiskCalculationForm.qml
    qml/features/risk/RiskSourcePanel.qml
    qml/features/risk/RiskRunView.qml
    qml/features/atmosphere/AtmosphereForm.qml
    qml/features/atmosphere/AtmosphereNumberField.qml
    qml/features/atmosphere/AtmosphereRangeField.qml
    qml/features/atmosphere/AtmosphereSourcePanel.qml
    qml/features/atmosphere/AtmosphereRunView.qml
    qml/controls/CustomSlider.qml
    qml/features/orbit/OrbitPage.qml
    qml/features/orbit/OrbitScene.qml
    qml/features/orbit/OrbitGrid.qml
    qml/features/orbit/OrbitLegend.qml
    qml/features/orbit/OrbitPath.qml
    qml/features/orbit/OrbitBody.qml
    qml/features/orbit/OrbitPlayback.qml
    qml/features/orbit/OrbitLayerDialog.qml
    qml/controls/CustomSwitch.qml
    qml/controls/CustomComboBox.qml
    qml/controls/CustomNumberField.qml
    qml/features/settings/CalculationSettingsDialog.qml
    qml/features/settings/SettingsPane.qml
    qml/features/settings/SettingRow.qml
    qml/features/settings/SettingsAppearancePane.qml
    qml/features/settings/ThemeTile.qml
    qml/features/settings/CalculationFitPane.qml
    qml/features/settings/CalculationPhotometryPane.qml
    qml/features/settings/CalculationIntegratorPane.qml
    qml/features/settings/CalculationPerturbationsPane.qml
    qml/features/settings/CalculationNonGravPane.qml
    qml/features/settings/CalculationAboutPane.qml
    qml/features/settings/DataSourcesDialog.qml
    qml/controls/DialogBackdrop.qml
    qml/controls/CustomProgressBar.qml
    qml/components/shared/KeyValueList.qml
    qml/shell/CalculationDialog.qml
    qml/shell/ImportObservationsDialog.qml
    qml/shell/RecentObservationFiles.qml
    qml/features/observations/ObservationDetails.qml
    qml/features/observations/ObservationDetailsDialog.qml
    qml/controls/InteractionCursor.qml
    qml/controls/CustomScrollBar.qml
    qml/controls/CustomIcon.qml
    qml/controls/CustomTab.qml
    qml/controls/CustomSegmentedControl.qml
    qml/controls/CustomTextLink.qml
    qml/components/shared/PlotLegendItem.qml
    qml/components/shared/InlineValue.qml
    qml/components/shared/ParameterTable.qml
    qml/shell/NavigationItem.qml
    qml/shell/SolutionHeader.qml
    qml/features/solution/ResidualSection.qml
    qml/components/shared/DiagnosticAxis.qml
    qml/features/solution/ResidualTimeChart.qml
    qml/features/solution/ResidualHistogram.qml
    qml/features/solution/SelectedObservationStrip.qml
    qml/components/shared/ResidualOffset.qml
    qml/features/solution/StationDiagnostics.qml
    qml/features/solution/StationResidualTable.qml
    qml/features/solution/StationResidualRow.qml
    qml/features/solution/OrbitDetailsDialog.qml
    qml/features/solution/MoidDetailsDialog.qml
    qml/features/solution/MoidSummary.qml
    qml/features/solution/DerivedOrbitParameters.qml
    qml/features/solution/CartesianStateTable.qml
    qml/features/solution/EstimatedSize.qml
    qml/shell/ObservationDropArea.qml
    qml/shared/UI.js
    qml/shared/AppTheme.qml
    qml/shared/DialogState.qml
    qml/controls/CustomText.qml
    qml/controls/CustomButton.qml
    qml/controls/CustomCheckBox.qml
    qml/controls/CustomTextField.qml
    qml/controls/BasePage.qml
    qml/controls/BaseDialog.qml
    qml/components/shared/SectionHeading.qml
    qml/components/shared/MetricValue.qml
    qml/shell/WorkspaceWindow.qml
    qml/shell/NavigationRail.qml
    qml/shell/CommandBar.qml
    qml/features/solution/SolutionPage.qml
    qml/features/solution/SolutionCharts.qml
    qml/features/solution/SolutionSummary.qml
    qml/features/solution/OrbitElementsTable.qml
    qml/features/solution/NonGravParameters.qml
    qml/features/solution/ResidualPlot.qml
    qml/features/solution/ResidualPointHint.qml
    qml/features/solution/ResidualReferenceCircle.qml
    qml/features/solution/CovariancePanel.qml
    qml/features/solution/CovarianceMatrix.qml
    qml/features/solution/CovarianceDetails.qml
    qml/components/shared/CustomMatrixCell.qml
    qml/features/observations/ObservationsPage.qml
    qml/features/observations/ObservationFilters.qml
    qml/features/observations/ObservationReviewDialog.qml
    qml/features/observations/ObservationHeaderCell.qml
    qml/features/observations/ObservationDetailSection.qml
    qml/features/observations/ObservationsTable.qml
    qml/features/observations/ObservationCell.qml
    qml/features/observations/ObservationInspector.qml
    qml/features/ephemerides/EphemeridesPage.qml
    qml/features/ephemerides/EphemerisSettingsDialog.qml
)

set(CLIENT_ICON_FILES
    resources/icons/workspace/search.svg
    resources/icons/workspace/plus.svg
    resources/icons/workspace/minus.svg
    resources/icons/workspace/rotate-ccw.svg
    resources/icons/workspace/view.svg
    resources/icons/workspace/pause.svg
    resources/icons/workspace/panel-left-close.svg
    resources/icons/workspace/panel-left-open.svg
    resources/app.ico
    resources/icons/workspace/app-logo.svg
    resources/icons/workspace/sliders-horizontal.svg
    resources/icons/workspace/activity.svg
    resources/icons/workspace/cloud-sun.svg
    resources/icons/workspace/sun.svg
    resources/icons/workspace/sigma.svg
    resources/icons/workspace/info.svg
    resources/icons/workspace/triangle-alert.svg
    resources/icons/workspace/check.svg
    resources/icons/workspace/orbit.svg
    resources/icons/workspace/scan-line.svg
    resources/icons/workspace/table-2.svg
    resources/icons/workspace/folder-open.svg
    resources/icons/workspace/folder.svg
    resources/icons/workspace/file-text.svg
    resources/icons/workspace/play.svg
    resources/icons/workspace/calendar-range.svg
    resources/icons/workspace/download.svg
    resources/icons/workspace/square.svg
    resources/icons/workspace/chevron-right.svg
    resources/icons/workspace/circle-dot.svg
    resources/icons/workspace/arrow-up-right.svg
    resources/icons/workspace/x.svg
    resources/icons/workspace/moon.svg
    resources/icons/workspace/LICENSE
)

set_source_files_properties(main.cpp src/features/settings/CalculationSettingsModel.cpp PROPERTIES
    COMPILE_DEFINITIONS "TRAJECTA_VERSION=\"${TRAJECTA_VERSION}\"")

if(WIN32)
    set(CLIENT_FONT_FILES
        resources/fonts/Inter.ttc
        resources/fonts/JetBrainsMono-Regular.ttf
        resources/fonts/JetBrainsMono-Bold.ttf)
else()
    set(CLIENT_FONT_FILES
        resources/fonts/SourceSans3-Regular.ttf
        resources/fonts/SourceSans3-It.ttf
        resources/fonts/SourceSans3-Semibold.ttf
        resources/fonts/SourceSans3-Bold.ttf
        resources/fonts/CascadiaMono.ttf)
endif()

qt_add_qml_module(TrajectaUI
    URI TrajectaUI
    VERSION 1.0
    SOURCES ${CLIENT_SOURCES}
    QML_FILES ${CLIENT_QML_FILES}
    RESOURCES
        ${CLIENT_ICON_FILES}
        ${CLIENT_FONT_FILES}
)

target_include_directories(TrajectaUI PRIVATE
    src
    src/features/orbit
    src/app
    src/features/solution
    src/features/observations
    src/controllers
    src/geometry
    src/io
    src/models
    src/settings
    src/util
)

target_link_libraries(TrajectaUI PRIVATE
    Qt6::Quick
    Qt6::QuickControls2
    Qt6::Quick3D
    Qt6::Charts
    Qt6::Widgets
    Qt6::Concurrent
    Qt6::Network
    OrbitalCore
    StrewnCore
)

# Updating the DLL need not relink the UI when its import library is unchanged.
if(WIN32)
    target_sources(TrajectaUI PRIVATE src/app/TaskbarProgress.h src/app/TaskbarProgress.cpp)
    target_link_libraries(TrajectaUI PRIVATE ole32 uuid)
    find_package(GDAL CONFIG REQUIRED)
    add_custom_target(TrajectaUI-deploy-core
        COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:TrajectaUI>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:OrbitalCore>
            $<TARGET_FILE_DIR:TrajectaUI>
        DEPENDS OrbitalCore
        COMMAND_EXPAND_LISTS
        VERBATIM
    )
    add_dependencies(TrajectaUI TrajectaUI-deploy-core)
    add_custom_target(TrajectaUI-deploy-strewn
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:StrewnCore> $<TARGET_FILE_DIR:TrajectaUI>
        COMMAND ${CMAKE_COMMAND}
            "-DRASTER_DLL=$<TARGET_FILE:GDAL::GDAL>"
            "-DDESTINATION=$<TARGET_FILE_DIR:TrajectaUI>"
            "-DCMAKE_GET_RUNTIME_DEPENDENCIES_COMMAND=${CMAKE_DUMPBIN}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/../StrewnField/cmake/DeployRaster.cmake"
        DEPENDS StrewnCore
        VERBATIM)
    add_dependencies(TrajectaUI TrajectaUI-deploy-strewn)
endif()

# DE440 kernels, observatory list, debias table — deployed from OrbitalCore/data
orbital_core_deploy_data(TrajectaUI)

# --- Tests ---
find_package(Qt6 REQUIRED COMPONENTS Test Qml)
enable_testing()

qt_add_executable(Tests
    ${CLIENT_SOURCES}
    tests/cpp/tst_workspace.h
    tests/cpp/tst_main.cpp
    tests/cpp/tst_recentfiles.h
    tests/cpp/tst_orbitgeometry.h
    tests/cpp/tst_cometorbits.h
    tests/cpp/tst_gridgeometry.h
    tests/cpp/tst_diagnosticsmodel.h
    tests/cpp/tst_nongravsettings.h
    tests/cpp/tst_settingsstorage.h
    tests/cpp/tst_orbitengine.h
)

target_include_directories(Tests PRIVATE
    src
    src/features/orbit
    src/app
    src/features/solution
    src/features/observations
    src/controllers
    src/geometry
    src/io
    src/models
    src/settings
    src/util
)
target_link_libraries(Tests PRIVATE Qt6::Test Qt6::Core Qt6::Qml Qt6::Gui Qt6::QuickControls2 Qt6::Quick3D Qt6::Charts Qt6::Widgets Qt6::Concurrent Qt6::Network OrbitalCore StrewnCore)
add_test(NAME Tests COMMAND Tests -o -,txt)
target_compile_definitions(Tests PRIVATE ORBITAL_UI_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
qt_add_resources(Tests "client_test_assets"
    PREFIX "/qt/qml/TrajectaUI"
    FILES ${CLIENT_FONT_FILES} ${CLIENT_ICON_FILES})
