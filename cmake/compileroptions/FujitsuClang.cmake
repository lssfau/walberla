message(STATUS "Setting FujitsuClang specific compiler options")

add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow -Wno-c++11-extensions -Qunused-arguments" )

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )

if ( WALBERLA_BUILD_WITH_LIKWID_MARKERS )
    find_library( LIKWID_LIB likwid HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )
    find_path( LIKWID_INCLUDE_DIR likwid.h HINTS $ENV{LIKWID_INCDIR} $ENV{LIKWID_ROOT}/include )

    # For some reason, these turned out to be necessary when building with likwid on Fugaku
     find_library( LIKWIDLUA_LIB NAMES likwid-lua HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )
     find_library( LIKWIDHWLOC_LIB NAMES likwid-hwloc HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )

    if ( LIKWID_LIB AND LIKWID_INCLUDE_DIR)
        set( LIKWID_FOUND 1 )
        include_directories( ${LIKWID_INCLUDE_DIR})
        add_definitions ( "-DLIKWID_PERFMON" )
        list ( APPEND SERVICE_LIBS ${LIKWID_LIB} )
            list ( APPEND SERVICE_LIBS ${LIKWIDLUA_LIB} ${LIKWIDHWLOC_LIB} )
    else()
        message(WARNING "likwid marker library not found. Set environment variable LIKWID_ROOT")
        set ( WALBERLA_BUILD_WITH_LIKWID_MARKERS OFF CACHE BOOL "Compile in markers for likwid-perfctr" FORCE )
    endif()
endif()