message(STATUS "Setting Cray specific compiler options")

# Fixes linker errors with Cray compiler
add_flag(CMAKE_EXE_LINKER_FLAGS "-dynamic -L/opt/gcc/4.9.3/snos/lib64")

# Silences compiler and linker warnings and information with the Cray compiler
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
add_flag(CMAKE_CXX_FLAGS "-h nomessage=1") # CC-1    The source file does not
                                           # end with a new-line character.
add_flag(CMAKE_C_FLAGS "-DSQLITE_HAVE_ISNAN") # SQLite will not work correctly
                                              # with the -ffast-math option of
                                              # GCC.
add_flag(CMAKE_CXX_FLAGS "-DSQLITE_HAVE_ISNAN") # SQLite will not work correctly
                                                # with the -ffast-math option of
                                                # GCC.

if(NOT WALBERLA_BUILD_WITH_OPENMP)
  add_flag(CMAKE_C_FLAGS "-h noomp")
  add_flag(CMAKE_CXX_FLAGS "-h noomp")
endif()

# Treat warnings as errors leaving it since it is not supported by
# COMPILE_WARNING_AS_ERROR (cmake 3.24)
if(WARNING_ERROR)
  add_flag(CMAKE_CXX_FLAGS "-h error_on_warning")
endif()
