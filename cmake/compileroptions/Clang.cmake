message(STATUS "Setting Clang specific compiler options")

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
   set(_apple_clang_minimal_version 11.0.0)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${_apple_clang_minimal_version})
    message(FATAL_ERROR "Clang version must be at least ${_apple_clang_minimal_version}!")
  endif()
else()
  set(_clang_minimal_version 7)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${_clang_minimal_version})
    message(FATAL_ERROR "Clang version must be at least ${_clang_minimal_version}!")
  endif()
endif()

if(WALBERLA_OPTIMIZE_FOR_LOCALHOST)

  if((CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID
                                                     STREQUAL "Clang")
     AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
    # no -march=native available on this compiler, but there is currently only
    # one such processor
  else()
    add_flag(CMAKE_CXX_FLAGS "-march=native")
    add_flag(CMAKE_C_FLAGS "-march=native")
  endif()

  if(EXISTS "/proc/sys/abi/sve_default_vector_length")
    file(READ "/proc/sys/abi/sve_default_vector_length" SVE_LENGTH_BYTES)
    string(STRIP "${SVE_LENGTH_BYTES}" SVE_LENGTH_BYTES)
    math(EXPR SVE_LENGTH "${SVE_LENGTH_BYTES} * 8")
    add_flag(CMAKE_CXX_FLAGS "-msve-vector-bits=${SVE_LENGTH}")
    add_flag(CMAKE_C_FLAGS "-msve-vector-bits=${SVE_LENGTH}")
  endif()
endif()

if(NOT WARNING_DEPRECATED)
  add_flag(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
endif()

add_flag(CMAKE_CXX_FLAGS
         "-Wall -Wconversion -Wshadow -Wno-c++11-extensions -Qunused-arguments")

if(WALBERLA_STL_BOUNDS_CHECKS)
  add_definitions("-D_GLIBCXX_DEBUG")
  add_definitions("-D_LIBCPP_DEBUG=1")
endif()

if(WALBERLA_BUILD_WITH_FASTMATH)
  add_flag(CMAKE_CXX_FLAGS "-ffast-math")
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
endif()

if(WALBERLA_SANITIZE_ADDRESS)
  add_flag(CMAKE_CXX_FLAGS "-fsanitize=address")
endif()

if(WALBERLA_SANITIZE_UNDEFINED)
  add_flag(CMAKE_CXX_FLAGS "-fsanitize=undefined")
endif()

if(WALBERLA_BUILD_WITH_OPENMP)
  if(APPLE
     AND EXISTS /opt/local/lib/libomp
     AND EXISTS /opt/local/include/libomp) # find libomp from MacPorts
    set(CMAKE_FRAMEWORK_PATH /opt/local/lib/libomp)
    set(CMAKE_INCLUDE_PATH /opt/local/include/libomp)
  endif()
  find_package(OpenMP)
  if(OpenMP_FOUND)
    add_flag(CMAKE_C_FLAGS "${OpenMP_C_FLAGS}")
    add_flag(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
    list(APPEND CORE_SERVICE_LIBS ${OpenMP_CXX_LIBRARIES})
    if(OpenMP_CXX_INCLUDE_DIRS)
      include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    endif()
  else()
    message(FATAL_ERROR "Could NOT enable OpenMP")
  endif()

  # check for bug in combination with OpenMP and sign conversion
  # https://bugs.llvm.org/show_bug.cgi?id=48387
  try_compile(
    WALBERLA_CLANG_OPENMP_BUG "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestClangOpenMPBug.cpp"
    COMPILE_DEFINITIONS -Werror)
  if(NOT ${WALBERLA_CLANG_OPENMP_BUG})
    message(
      WARNING
        "Setting -Wno-sign-conversion due to a compiler bug in LLVM (https://bugs.llvm.org/show_bug.cgi?id=48387)"
    )
    add_flag(CMAKE_CXX_FLAGS "-Wno-sign-conversion")
  endif()
endif()