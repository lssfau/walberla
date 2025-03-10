message(STATUS "Setting IBM specific compiler options")

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.1.1)
  message(FATAL_ERROR "IBM compiler version must be at least 17.1.1!")
endif()

# Fixes linker errors with IBM compiler
add_flag(CMAKE_CXX_FLAGS "-qpic=large")

# Silences compiler and linker warnings and information with the IBM compiler
add_flag(CMAKE_CXX_FLAGS "-qsuppress=1586-267") # 1586-267 (I) Inlining of
                                                # specified subprogram failed
                                                # due to the presence of a C++
                                                # exception handler
add_flag(CMAKE_CXX_FLAGS "-qsuppress=1586-266") # 1586-266 (I) Inlining of
                                                # specified subprogram failed
                                                # due to the presence of a
                                                # global label
add_flag(CMAKE_CXX_FLAGS "-qsuppress=1500-030") # 1500-030: (I) INFORMATION:
                                                # [...] Additional optimization
                                                # may be attained by recompiling
                                                # and specifying MAXMEM option
                                                # with a value greater than
                                                # 8192.
add_flag(CMAKE_C_FLAGS "-qsuppress=1500-030") # 1500-030: (I) INFORMATION: [...]
                                              # Additional optimization may be
                                              # attained by recompiling and
                                              # specifying MAXMEM option with a
                                              # value greater than 8192.
