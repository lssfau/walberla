message(STATUS "Setting IntelLLVM specific compiler options")

# fastmath
if(NOT WALBERLA_BUILD_WITH_FASTMATH)
  add_flag(CMAKE_CXX_FLAGS "-fp-model=precise")
endif()

set(CMAKE_C_FLAGS_DEBUGOPTIMIZED "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3")
set(CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3")

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
    list(APPEND SERVICE_LIBS ${OpenMP_CXX_LIBRARIES})
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
