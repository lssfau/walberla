message(STATUS "Setting IntelLLVM specific compiler options")

#fastmath
if ( NOT WALBERLA_BUILD_WITH_FASTMATH )
   add_flag( CMAKE_CXX_FLAGS "-fp-model=precise")
endif()

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )