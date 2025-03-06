message(STATUS "Setting GNU specific compiler options")

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
   message(FATAL_ERROR "GCC version must be at least 8!")
endif()

# Profile guided optimization
if ( WALBERLA_PROFILE_GENERATE )
   add_flag( CMAKE_CXX_FLAGS "-fprofile-generate" )
endif()
if ( WALBERLA_PROFILE_USE )
   add_flag( CMAKE_CXX_FLAGS "-fprofile-use" )
endif()


# common flags for g++
if ( NOT WARNING_DISABLE )
   add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow" )
endif()

# architecture optimization
if( WALBERLA_OPTIMIZE_FOR_LOCALHOST )
   add_flag ( CMAKE_CXX_FLAGS "-march=native" )
endif()

if( EXISTS "/proc/sys/abi/sve_default_vector_length" )
   file( READ "/proc/sys/abi/sve_default_vector_length" SVE_LENGTH_BYTES )
   string(STRIP "${SVE_LENGTH_BYTES}" SVE_LENGTH_BYTES)
   math(EXPR SVE_LENGTH "${SVE_LENGTH_BYTES} * 8")
   add_flag ( CMAKE_CXX_FLAGS "-msve-vector-bits=${SVE_LENGTH}" )
   add_flag ( CMAKE_C_FLAGS   "-msve-vector-bits=${SVE_LENGTH}" )
endif()

# Warning flags
add_flag ( CMAKE_CXX_FLAGS "-Wfloat-equal -Wextra" )

if ( WARNING_PEDANTIC )
   add_flag ( CMAKE_CXX_FLAGS "-pedantic" )
endif ( )

 # omit deprecated warnings
if( NOT WARNING_DEPRECATED)
       add_flag ( CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
   endif()

if ( WALBERLA_STL_BOUNDS_CHECKS )
   add_definitions ( "-D_GLIBCXX_DEBUG" )
   add_definitions ( "-D_LIBCPP_DEBUG=1" )
endif()

if( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0" )
   # Omit maybe-uninitialized for gcc 12 for now. Check if it is a bug or a real problem:
   add_flag( CMAKE_CXX_FLAGS "-Wno-maybe-uninitialized" )
   # GCC 12 reports a "array bounds" warning at UniformBufferedScheme.h:297 (error: array subscript 26 is above array bounds of)
   # when e.g. compiling the GhostLayerCommTest.
   # Since this is most probably a bug in GCC disable the warning for now
   add_flag( CMAKE_CXX_FLAGS "-Wno-array-bounds" )
endif()

if ( WALBERLA_BUILD_WITH_FASTMATH )
   add_flag( CMAKE_CXX_FLAGS "-ffast-math")
endif()

#GCC 5+ ABI selection
if( NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0.0 )
   option ( WALBERLA_USE_CPP11_ABI "On GCC 5+ use the C++11 ABI" ON )
   if( WALBERLA_USE_CPP11_ABI )
      add_flag( CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=1" )
   else()
      add_flag( CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0" )
   endif()
endif()

if (NOT WIN32)
find_package ( Backtrace QUIET )
if ( Backtrace_FOUND )
   list ( APPEND SERVICE_LIBS ${Backtrace_LIBRARIES} )
   set ( WALBERLA_BUILD_WITH_BACKTRACE ON )
   set ( WALBERLA_BACKTRACE_HEADER ${Backtrace_HEADER} )
endif ( Backtrace_FOUND )
endif()

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )

# old nvcc compilers and newer stdlibc++ are incompatible. This needs to be checked!
if (WALBERLA_STL_BOUNDS_CHECKS AND WALBERLA_BUILD_WITH_CODEGEN AND WALBERLA_BUILD_WITH_CUDA AND CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "11.0" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
   message(FATAL_ERROR "WALBERLA_STL_BOUNDS_CHECKS is not compatible with your CUDA compiler")
endif()

if ( WALBERLA_BUILD_WITH_GPROF )
   add_flag ( CMAKE_CXX_FLAGS        "-pg" )
endif()

if ( WALBERLA_SANITIZE_ADDRESS )
   add_flag( CMAKE_CXX_FLAGS "-fsanitize=address")
endif()

if ( WALBERLA_SANITIZE_UNDEFINED )
   add_flag( CMAKE_CXX_FLAGS "-fsanitize=undefined")
endif()

############################################################################################################################
##
##  Testing Coverage
##
############################################################################################################################
if (WALBERLA_BUILD_WITH_GCOV AND CMAKE_COMPILER_IS_GNUCXX  )
    add_flag ( CMAKE_CXX_FLAGS "--coverage" )
endif()
############################################################################################################################