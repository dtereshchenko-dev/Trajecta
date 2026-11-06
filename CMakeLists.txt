cmake_minimum_required(VERSION 3.21)
include(cmake/TrajectaVersion.cmake)
project(Trajecta VERSION ${TRAJECTA_VERSION} LANGUAGES C CXX)

if(MSVC AND CMAKE_GENERATOR MATCHES "Ninja")
    set(_msvc_probe "${CMAKE_BINARY_DIR}/CMakeFiles/OrbitalShowIncludes")
    file(MAKE_DIRECTORY "${_msvc_probe}")
    set(_msvc_environment "${_msvc_probe}/compile.cmd")
    file(WRITE "${_msvc_environment}" "@echo off\r\nchcp 65001 >nul\r\nset VSLANG=1033\r\n%*\r\nexit /b %errorlevel%\r\n")
    file(WRITE "${_msvc_probe}/probe.h" "\n")
    file(WRITE "${_msvc_probe}/probe.c" "#include \"probe.h\"\nint probe;\n")
    execute_process(COMMAND ${_msvc_environment} "${CMAKE_C_COMPILER}" /nologo /showIncludes /c probe.c
        WORKING_DIRECTORY "${_msvc_probe}"
        OUTPUT_VARIABLE _msvc_output ERROR_VARIABLE _msvc_error
        RESULT_VARIABLE _msvc_result ENCODING NONE)
    if(NOT _msvc_result EQUAL 0 OR NOT _msvc_output MATCHES
            "(^|\n)([^:\n][^:\n]+:[^:\n]*[^: \n][^: \n]:?[ \t]+)([A-Za-z]:\\\\|/)")
        message(FATAL_ERROR "Cannot determine MSVC header dependency prefix: ${_msvc_output}${_msvc_error}")
    endif()
    set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_MATCH_2}")
    set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "${CMAKE_MATCH_2}")
    set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "${CMAKE_MATCH_2}")
    list(PREPEND CMAKE_C_COMPILER_LAUNCHER ${_msvc_environment})
    list(PREPEND CMAKE_CXX_COMPILER_LAUNCHER ${_msvc_environment})
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    foreach(cfg RELEASE RELWITHDEBINFO MINSIZEREL)
        string(REGEX REPLACE "/Ob[0-9]" "" CMAKE_CXX_FLAGS_${cfg}
               "${CMAKE_CXX_FLAGS_${cfg}}")
        string(STRIP "${CMAKE_CXX_FLAGS_${cfg}}" CMAKE_CXX_FLAGS_${cfg})
    endforeach()
endif()

enable_testing()

add_subdirectory(Common/Http)
add_subdirectory(OrbitalCore)
add_subdirectory(TrajectaCLI)
add_subdirectory(StrewnField)
add_subdirectory(TrajectaUI)
