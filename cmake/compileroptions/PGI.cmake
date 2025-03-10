message(STATUS "Setting PGI specific compiler options")

# Silences compiler and linker warnings and information with the PGI compiler
add_flag(CMAKE_CXX_FLAGS "--display_error_number")
add_flag(CMAKE_C_FLAGS "--display_error_number")
if(CMAKE_VERSION VERSION_LESS "3.19")
  # https://github.com/Kitware/CMake/commit/52eee1938919deb59cc2b51d44f365f0d9a418e5
  set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION
      "--c++${CMAKE_CXX_STANDARD}")
endif()
add_flag(CMAKE_CXX_FLAGS "--diag_suppress=1") # last line of file ends without a
                                              # newline
add_flag(CMAKE_CXX_FLAGS "--diag_suppress=111") # statement is unreachable
add_flag(CMAKE_C_FLAGS "--diag_suppress=111") # statement is unreachable
add_flag(CMAKE_C_FLAGS "--diag_suppress=550") # variable [...] was set but never
                                              # used
add_flag(CMAKE_C_FLAGS "--diag_suppress=191") # type qualifier is meaningless on
                                              # cast type
