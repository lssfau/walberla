message(STATUS "Setting Intel specific compiler options")

set(_intel_minimal_version 19.0.0)
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${_intel_minimal_version})
  message(FATAL_ERROR "IBM compiler version must be at least ${_intel_minimal_version}!")
endif()

if(WALBERLA_PROFILE_GENERATE)
  add_flag(CMAKE_CXX_FLAGS "-prof-gen")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profile")
  add_flag(CMAKE_CXX_FLAGS "-prof-dir${CMAKE_BINARY_DIR}/profile")
endif()

if(WALBERLA_PROFILE_USE)
  add_flag(CMAKE_CXX_FLAGS "-prof-use")
  add_flag(CMAKE_CXX_FLAGS "-prof-dir${CMAKE_BINARY_DIR}/profile")
endif()

# common flags for intel and g++
if(NOT WARNING_DISABLE)
  add_flag(CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow")
endif()

# architecture optimization
if(WALBERLA_OPTIMIZE_FOR_LOCALHOST)
  add_flag(CMAKE_CXX_FLAGS "-march=native")
  add_flag(CMAKE_C_FLAGS "-march=native")

  add_flag(CMAKE_CXX_FLAGS "-xhost")
  add_flag(CMAKE_C_FLAGS "-xhost")

  if(EXISTS "/proc/sys/abi/sve_default_vector_length")
    file(READ "/proc/sys/abi/sve_default_vector_length" SVE_LENGTH_BYTES)
    string(STRIP "${SVE_LENGTH_BYTES}" SVE_LENGTH_BYTES)
    math(EXPR SVE_LENGTH "${SVE_LENGTH_BYTES} * 8")
    add_flag(CMAKE_CXX_FLAGS "-msve-vector-bits=${SVE_LENGTH}")
    add_flag(CMAKE_C_FLAGS "-msve-vector-bits=${SVE_LENGTH}")
  endif()
endif()

# system headers are also supported by intel, but cmake does not recognize that
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
add_flag(CMAKE_CXX_FLAGS "-wd2928,2504,2259,1682,597")
# disable icc/icpc deprecation warning
add_flag(CMAKE_CXX_FLAGS "-diag-disable=10441")

# omit deprecated warnings
if(NOT WARNING_DEPRECATED)
  add_flag(CMAKE_CXX_FLAGS "-wd1478") # Disable compiler warning # 1478:
                                      # "declared as deprecated"
endif()

if(WALBERLA_STL_BOUNDS_CHECKS)
  add_definitions("-D_GLIBCXX_DEBUG")
  add_definitions("-D_LIBCPP_DEBUG=1")
endif()

if(WALBERLA_BUILD_WITH_FASTMATH)
  add_flag(CMAKE_CXX_FLAGS "-fp-model fast=2 -no-prec-sqrt -no-prec-div")
endif()

find_package(Backtrace QUIET)
if(Backtrace_FOUND)
 list(APPEND CORE_SERVICE_LIBS ${Backtrace_LIBRARIES})
 set(WALBERLA_BUILD_WITH_BACKTRACE ON)
 set(WALBERLA_BACKTRACE_HEADER ${Backtrace_HEADER})
endif(Backtrace_FOUND)

set(CMAKE_C_FLAGS_DEBUGOPTIMIZED "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3")
set(CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3")

if(WALBERLA_BUILD_WITH_GPROF)
  add_flag(CMAKE_CXX_FLAGS "-pg")
  add_flag(CMAKE_EXE_LINKER_FLAGS "-pg")
endif()