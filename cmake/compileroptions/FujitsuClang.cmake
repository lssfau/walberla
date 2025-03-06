message(STATUS "Setting FujitsuClang specific compiler options")

add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow -Wno-c++11-extensions -Qunused-arguments" )

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )

if ( WALBERLA_BUILD_WITH_LIKWID_MARKERS )

   # For some reason, these turned out to be necessary when building with likwid on Fugaku
   find_library( LIKWIDLUA_LIB NAMES likwid-lua HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )
   find_library( LIKWIDHWLOC_LIB NAMES likwid-hwloc HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )

   if ( LIKWID_LIB AND LIKWID_INCLUDE_DIR)
      list ( APPEND SERVICE_LIBS ${LIKWIDLUA_LIB} ${LIKWIDHWLOC_LIB} )
   endif()
endif()