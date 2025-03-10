message(STATUS "Setting NEC specific compiler options")

# C++ language features for NEC compiler
set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION
    "-Kcpp${CMAKE_CXX_STANDARD}")
set(CMAKE_TRY_COMPILE_PLATFORM_VARIABLES
    CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION)
add_flag(
  CMAKE_CXX_FLAGS
  "${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} -Krtti -Kexceptions -size_t64 -Kgcc"
)
add_flag(CMAKE_CXX_FLAGS "-D__BIG_ENDIAN -D__BYTE_ORDER=__BIG_ENDIAN")
add_flag(CMAKE_CXX_FLAGS "-Tnoauto,used")
add_flag(CMAKE_EXE_LINKER_FLAGS "-Wl,-h,muldefs")
add_flag(CMAKE_C_FLAGS "-size_t64 -Kgcc")
add_flag(CMAKE_C_FLAGS "-D__BIG_ENDIAN -D__BYTE_ORDER=__BIG_ENDIAN")
add_flag(CMAKE_C_FLAGS "-DSQLITE_OMIT_WAL -DHAVE_UTIME -DTHREADSAFE=0")
set(CMAKE_RANLIB /bin/true)
set(CMAKE_SKIP_BUILD_RPATH TRUE)
set(CMAKE_C_FLAGS_DEBUGOPTIMIZED "-Chopt -g")
set(CMAKE_C_FLAGS_DEBUG "-Cdebug -g")
set(CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "-Chopt -g")
set(CMAKE_CXX_FLAGS_DEBUG "-Cdebug -g")

add_flag(CMAKE_CXX_FLAGS "-wall")

# ##############################################################################
#
# Fix compiler bugs
#
# ##############################################################################

# The NEC SX has a few issues in its standard library headers
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/math.h
  "#include_next <math.h>\n#undef fpclassify\n#undef signbit\n#undef isfinite\n#undef isinf\n#undef isnan\n#undef isnormal\n#undef isgreater\n#undef isgreaterequal\n#undef isless\n#undef islessequal\n#undef islessgreater\n#undef isunordered\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/types.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/types.h\"   \n#undef uint_t\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/acl.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/acl.h\"     \n#undef uint_t\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/if_ehcpl.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/if_ehcpl.h\"\n#undef uint_t\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/ptms.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/ptms.h\"    \n#undef uint_t\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/stream.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/stream.h\"  \n#undef uint_t\n"
)
file(
  WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/strsubr.h
  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/strsubr.h\" \n#undef uint_t\n"
)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/math.h
               ${walberla_BINARY_DIR}/src/math.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/types.h
               ${walberla_BINARY_DIR}/src/sys/types.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/acl.h
               ${walberla_BINARY_DIR}/src/sys/acl.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/if_ehcpl.h
               ${walberla_BINARY_DIR}/src/sys/if_ehcpl.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/ptms.h
               ${walberla_BINARY_DIR}/src/sys/ptms.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/stream.h
               ${walberla_BINARY_DIR}/src/sys/stream.h COPYONLY)
configure_file(${walberla_BINARY_DIR}/CMakeFiles/src/sys/strsubr.h
               ${walberla_BINARY_DIR}/src/sys/strsubr.h COPYONLY)

if(WALBERLA_BUILD_WITH_OPENMP)
  message(STATUS "Enabling OpenMP workaround for NEC")
  add_flag(CMAKE_C_FLAGS "-Popenmp")
  add_flag(CMAKE_CXX_FLAGS "-Popenmp")
endif()
