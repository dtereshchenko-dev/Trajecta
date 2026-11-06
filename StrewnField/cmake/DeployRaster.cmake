get_filename_component(source_directory "${RASTER_DLL}" DIRECTORY)
set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "windows+pe")
set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "dumpbin")
file(GET_RUNTIME_DEPENDENCIES
    LIBRARIES "${RASTER_DLL}"
    DIRECTORIES "${source_directory}"
    RESOLVED_DEPENDENCIES_VAR dependencies
    UNRESOLVED_DEPENDENCIES_VAR unresolved
    PRE_EXCLUDE_REGEXES "[Aa][Pp][Ii]-[Mm][Ss]-" "[Ee][Xx][Tt]-[Mm][Ss]-"
    POST_EXCLUDE_REGEXES ".*[/\\][Ss][Yy][Ss][Tt][Ee][Mm]32[/\\].*")
if(unresolved)
    message(FATAL_ERROR "Unresolved GDAL runtime dependencies: ${unresolved}")
endif()
list(APPEND dependencies "${RASTER_DLL}")
foreach(dependency IN LISTS dependencies)
    get_filename_component(name "${dependency}" NAME)
    if(NOT "${dependency}" STREQUAL "${DESTINATION}/${name}")
        file(COPY_FILE "${dependency}" "${DESTINATION}/${name}" ONLY_IF_DIFFERENT)
    endif()
endforeach()
get_filename_component(source_prefix "${source_directory}" DIRECTORY)
if(NOT EXISTS "${source_prefix}/share/proj/proj.db" AND EXISTS "${source_prefix}/../share/proj/proj.db")
    get_filename_component(source_prefix "${source_prefix}" DIRECTORY)
endif()
if(EXISTS "${source_prefix}/share/proj/proj.db")
    file(MAKE_DIRECTORY "${DESTINATION}/share/proj")
    file(COPY_FILE "${source_prefix}/share/proj/proj.db" "${DESTINATION}/share/proj/proj.db" ONLY_IF_DIFFERENT)
endif()
if(EXISTS "${source_prefix}/share/gdal/grib2_table_4_2_0_0.csv")
    file(COPY "${source_prefix}/share/gdal/" DESTINATION "${DESTINATION}/share/gdal"
        PATTERN "*.cmake" EXCLUDE PATTERN "vcpkg*" EXCLUDE PATTERN "usage" EXCLUDE)
endif()
