cmake_minimum_required(VERSION 3.21)
project(TrajectaCLI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT TARGET OrbitalCore)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../OrbitalCore OrbitalCore_build)
endif()

add_executable(${PROJECT_NAME} main.cpp)
find_package(nlohmann_json CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OrbitalCore nlohmann_json::nlohmann_json)

if(WIN32)
    # Updating the DLL need not relink the CLI when its import library is unchanged.
    add_custom_target(${PROJECT_NAME}-deploy-core
        COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:OrbitalCore>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:OrbitalCore>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
        DEPENDS OrbitalCore
        COMMAND_EXPAND_LISTS
        VERBATIM
    )
    add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}-deploy-core)
endif()

orbital_core_deploy_data(${PROJECT_NAME})
