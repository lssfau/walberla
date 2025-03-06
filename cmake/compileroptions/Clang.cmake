message(STATUS "Setting Clang specific compiler options")
if( WALBERLA_OPTIMIZE_FOR_LOCALHOST )

   if(( CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64" )
      # no -march=native available on this compiler, but there is currently only one such processor
   else()
      add_flag ( CMAKE_CXX_FLAGS "-march=native" )
      add_flag ( CMAKE_C_FLAGS   "-march=native" )
   endif()

   if( EXISTS "/proc/sys/abi/sve_default_vector_length" )
      file( READ "/proc/sys/abi/sve_default_vector_length" SVE_LENGTH_BYTES )
      string(STRIP "${SVE_LENGTH_BYTES}" SVE_LENGTH_BYTES)
      math(EXPR SVE_LENGTH "${SVE_LENGTH_BYTES} * 8")
      add_flag ( CMAKE_CXX_FLAGS "-msve-vector-bits=${SVE_LENGTH}" )
      add_flag ( CMAKE_C_FLAGS   "-msve-vector-bits=${SVE_LENGTH}" )
   endif()
endif()

if( NOT WARNING_DEPRECATED)
   add_flag ( CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
endif()

add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow -Wno-c++11-extensions -Qunused-arguments" )

if ( WALBERLA_STL_BOUNDS_CHECKS )
   add_definitions ( "-D_GLIBCXX_DEBUG" )
   add_definitions ( "-D_LIBCPP_DEBUG=1" )
endif()

if ( WALBERLA_BUILD_WITH_FASTMATH )
   add_flag( CMAKE_CXX_FLAGS "-ffast-math")
endif()

if ( NOT WIN32 )
   find_package ( Backtrace QUIET )
   if ( Backtrace_FOUND )
      list ( APPEND SERVICE_LIBS ${Backtrace_LIBRARIES} )
      set ( WALBERLA_BUILD_WITH_BACKTRACE ON )
      set ( WALBERLA_BACKTRACE_HEADER ${Backtrace_HEADER} )
   endif ( Backtrace_FOUND )
endif()

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )

if ( WALBERLA_BUILD_WITH_GPROF )
   add_flag ( CMAKE_CXX_FLAGS        "-pg" )
endif()

if ( WALBERLA_SANITIZE_ADDRESS )
   add_flag( CMAKE_CXX_FLAGS "-fsanitize=address")
endif()

if ( WALBERLA_SANITIZE_UNDEFINED )
   add_flag( CMAKE_CXX_FLAGS "-fsanitize=undefined")
endif()