############################################################################################################################
##
## waLBerla's main cmake file
##
## Contents:
##   - definition of build options
##   - compiler variables ( c++ standard, warnings etc. )
##   - Finding of service libraries. Required: none, Optional: Boost, MPI, FFTW3, METIS, OpenMesh, Python
##     the include paths are set, and the libraries are added to variable SERVICE_LIBS
##   - Subdirectory cmake lists are called
##       -> src/   this folder contains all modules, each module (that contains c or cpp files) is linked to a
##                 static library.  Dependencies between these shared libraries are tracked manually,
##                 for more information see waLBerlaModuleDependencySystem.cmake
##       -> tests/ Same subdirectories as src/ folder. Contains tests for each module
##   - Export of variables into internal cache variables, for usage in applications or projects that use walberla as
##     subdirectory. Variables containing the service-libs,
##
############################################################################################################################



############################################################################################################################
##
## Project name, version, Custom CMake functions
##
############################################################################################################################

CMAKE_MINIMUM_REQUIRED (VERSION 3.10)


PROJECT ( walberla )

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${walberla_SOURCE_DIR}/cmake )

include ( waLBerlaFunctions )

set_version( 5 0 )

include( CMakeParseArguments )

# Enable CTest
enable_testing()
include( CTest )


############################################################################################################################




############################################################################################################################
##
## Definition of build options
##
############################################################################################################################


# Build options
option ( WALBERLA_DOUBLE_ACCURACY           "Floating point accuracy, defaults to double"     ON )
option ( WALBERLA_ENABLE_GUI                "Compile with GUI"                                   )

option ( WALBERLA_BUILD_TESTS               "Build Testcases"                                    )
option ( WALBERLA_BUILD_BENCHMARKS          "Build Benchmarks"                                ON )
option ( WALBERLA_BUILD_TOOLS               "Build Tools"                                        )
option ( WALBERLA_BUILD_TUTORIALS           "Build Tutorials"                                 ON )
option ( WALBERLA_BUILD_SHOWCASES           "Build Showcases"                                OFF )
option ( WALBERLA_BUILD_DOC                 "Build Documentation"                             ON )

option ( WALBERLA_BUILD_WITH_MPI            "Build with MPI"                                  ON )
option ( WALBERLA_BUILD_WITH_METIS          "Build with metis graph partitioner"             OFF )
option ( WALBERLA_BUILD_WITH_PARMETIS       "Build with ParMetis graph partitioner"          OFF )

option ( WALBERLA_BUILD_WITH_GPROF          "Enables gprof"                                      )
option ( WALBERLA_BUILD_WITH_GCOV           "Enables gcov"                                       )
option ( WALBERLA_BUILD_WITH_LTO            "Enable link time optimizations"                     )
option ( WALBERLA_BUILD_WITH_OPENMP         "Enable OpenMP support"                              )
option ( WALBERLA_BUILD_WITH_PYTHON         "Support for embedding Python"                       )
option ( WALBERLA_BUILD_WITH_CODEGEN        "Enable pystencils code generation"              OFF )


option ( WALBERLA_BUILD_WITH_LIKWID_MARKERS "Compile in markers for likwid-perfctr"              )

option ( WALBERLA_BUILD_WITH_CUDA	        "Enable CUDA support"                                )


option ( WALBERLA_BUILD_WITH_FASTMATH       "Fast math"                                          )

option ( WALBERLA_SIMD_FORCE_SCALAR         "Do not use SIMD operations even when available" OFF )

option ( WALBERLA_BUFFER_DEBUG              "Type checking for BufferSystem ( slow )"        OFF )

# Profile guided optimization
option ( WALBERLA_PROFILE_GENERATE          "Generates Profile for Optimization"                 )
option ( WALBERLA_PROFILE_USE               "Uses Profile to optimize"                           )

# Compiler Optimization
option ( WALBERLA_OPTIMIZE_FOR_LOCALHOST    "Enable compiler optimizations spcific to localhost" )

option ( WALBERLA_LOG_SKIPPED               "Log skipped cmake targets"                      ON  )

option ( WALBERLA_GIT_SUBMODULE_AUTO        "Check submodules during cmake run"               ON )

# Installation Directory
set ( CMAKE_INSTALL_PREFIX /usr/local/waLBerla CACHE STRING "The default installation directory."   )

# Default build type
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release CACHE STRING "Build Types: Debug Release DebugOptimized RelWithDebInfo MinSizeRel."  FORCE )
endif()
SET_PROPERTY( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release DebugOptimized RelWithDebInfo MinSizeRel )

# Debugging options                                      )
option ( WALBERLA_STL_BOUNDS_CHECKS  "Use debug capabilites of libstd++: iterator and bounds checks" )

# Warning options
option ( WARNING_DISABLE    "Disables additional compiler warnings"          OFF )
option ( WARNING_PEDANTIC   "Enables pedantic compiler warnings"             ON  )
option ( WARNING_ERROR      "Convert warnings to errors compiler warnings"   OFF )
option ( WARNING_DEPRECATED "Show warning when deprecated features are used" ON  )

# Sanitizer options
option ( WALBERLA_SANITIZE_ADDRESS    "Enables address sanitizer in gcc and clang"            )
option ( WALBERLA_SANITIZE_UNDEFINED  "Enables undefined behavior sanitizer in gcc and clang" )

# Every folder that is listed here can contain modules or tests
# this can be extended by applications to have own modules
# Here the src/ folder is added to this list, where all modules are located
list( APPEND WALBERLA_MODULE_DIRS "${walberla_SOURCE_DIR}/src" "${walberla_SOURCE_DIR}/tests" )
list( REMOVE_DUPLICATES  WALBERLA_MODULE_DIRS )
set ( WALBERLA_MODULE_DIRS  ${WALBERLA_MODULE_DIRS} CACHE INTERNAL "All folders that contain modules or tests" )

# target_link_libraries needs to called with keywords everywhere if it is called with keywords once
if( DEFINED CUDA_LINK_LIBRARIES_KEYWORD AND NOT CUDA_LINK_LIBRARIES_KEYWORD STREQUAL "" )
    set( WALBERLA_LINK_LIBRARIES_KEYWORD PUBLIC )
endif()

############################################################################################################################

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )

############################################################################################################################
##
## Compiler detection
##
############################################################################################################################

# Check for intel compiler
if( CMAKE_CXX_COMPILER MATCHES "icpc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "icpc" )
    option ( WALBERLA_CXX_COMPILER_IS_INTEL "Use Intel compiler" ON  )
    # Intel(R) Compiler has its own library archiver,
    # if you build libraries and do not use xiar,
    # the Intel compiler will complain about invalid
    # archives at the link phase.
    # The Intel(R) archiver is "xiar" usually
    # located in the same folder as the compiler,
    FIND_PROGRAM(XIAR xiar)
    IF(XIAR)
        SET(CMAKE_AR "${XIAR}")
    ENDIF(XIAR)
    MARK_AS_ADVANCED(XIAR)

    # Intel(R) Compiler also comes with its own linker
    # which provides a number of additional benefits when
    # linking code compiled with the Intel(R) compiler.
    # Again, usually in the same place as icc itself,
    FIND_PROGRAM(XILD xild)
    IF(XILD)
       SET(CMAKE_LINKER "${XILD}")
    ENDIF(XILD)
    MARK_AS_ADVANCED(XILD)

    if( CMAKE_VERSION VERSION_LESS 3.11.0 )
      set( CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17" )
      set( CMAKE_TRY_COMPILE_PLATFORM_VARIABLES CMAKE_CXX17_STANDARD_COMPILE_OPTION )
      add_flag ( CMAKE_CXX_FLAGS ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} )
    endif()
    if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.0.5" )
      # std::filesystem uses ABI tags, which don't work 19.0.2 but do in 19.0.5
      add_flag ( CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0" )
    endif()
else()
    option ( WALBERLA_CXX_COMPILER_IS_INTEL "Use Intel compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_INTEL )

# Check for Gnu compiler
if ( CMAKE_COMPILER_IS_GNUCXX  AND NOT WALBERLA_CXX_COMPILER_IS_INTEL )
     option ( WALBERLA_CXX_COMPILER_IS_GNU "Use gnu compiler" ON  )
else()
     option ( WALBERLA_CXX_COMPILER_IS_GNU "Use gnu compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_GNU )

# Check for Visual Studio
if ( MSVC )
     option ( WALBERLA_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" ON  )
else()
     option ( WALBERLA_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_MSVC )

# Check for IBM compiler
if( CMAKE_CXX_COMPILER MATCHES "xlc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "xlc" )
    option ( WALBERLA_CXX_COMPILER_IS_IBM "Use IBM compiler" ON  )
else()
    option ( WALBERLA_CXX_COMPILER_IS_IBM "Use IBM compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_IBM )

# Check for NEC SX compiler
if( CMAKE_CXX_COMPILER MATCHES "sxc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "sxc" OR CMAKE_CXX_COMPILER MATCHES "sxmpic" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "sxmpic" )
    option ( WALBERLA_CXX_COMPILER_IS_NEC "Use NEC compiler" ON  )
else()
    option ( WALBERLA_CXX_COMPILER_IS_NEC "Use NEC compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_NEC )

# Check for Clang compiler
if( CMAKE_CXX_COMPILER MATCHES "clang" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" )
    option ( WALBERLA_CXX_COMPILER_IS_CLANG "Use clang compiler" ON  )
else()
    option ( WALBERLA_CXX_COMPILER_IS_CLANG "Use clang compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_CLANG )

if( CMAKE_CXX_COMPILER_ID MATCHES Cray )
    option ( WALBERLA_CXX_COMPILER_IS_CRAY "Use Cray compiler" ON   )
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.4)
        message( FATAL_ERROR "Insufficient Cray Compiler Environment version" )
    endif()
else()
    option ( WALBERLA_CXX_COMPILER_IS_CRAY "Use Cray compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_CRAY )

# Check for MPI wrapper
get_filename_component( CXX_COMPILER_WITHOUT_PATH ${CMAKE_CXX_COMPILER} NAME )
if( CXX_COMPILER_WITHOUT_PATH MATCHES "mpi" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "mpi" )
    option ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" ON  )
else()
    option ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER )

############################################################################################################################



############################################################################################################################
##
## Compiler Setup
##
############################################################################################################################

# Profile guided optimization
if ( WALBERLA_PROFILE_GENERATE )
    if( WALBERLA_CXX_COMPILER_IS_INTEL )
        add_flag( CMAKE_CXX_FLAGS "-prof-gen" )
        file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profile" )
        add_flag( CMAKE_CXX_FLAGS "-prof-dir${CMAKE_BINARY_DIR}/profile" )
    elseif( WALBERLA_CXX_COMPILER_IS_GNU )
        add_flag( CMAKE_CXX_FLAGS "-fprofile-generate" )
    elseif( WALBERLA_CXX_COMPILER_IS_MSVC )
      add_flag ( CMAKE_CXX_FLAGS           "/GL"                )
      add_flag ( CMAKE_MODULE_LINKER_FLAGS "/LTCG:PGINSTRUMENT" )
      add_flag ( CMAKE_SHARED_LINKER_FLAGS "/LTCG:PGINSTRUMENT" )
      add_flag ( CMAKE_EXE_LINKER_FLAGS    "/LTCG:PGINSTRUMENT" )
    endif()
endif()

if ( WALBERLA_PROFILE_USE )
    if( WALBERLA_CXX_COMPILER_IS_INTEL )
       add_flag( CMAKE_CXX_FLAGS "-prof-use" )
       add_flag( CMAKE_CXX_FLAGS "-prof-dir${CMAKE_BINARY_DIR}/profile" )
    elseif( WALBERLA_CXX_COMPILER_IS_GNU )
       add_flag( CMAKE_CXX_FLAGS "-fprofile-use" )
    elseif( WALBERLA_CXX_COMPILER_IS_MSVC )
      add_flag ( CMAKE_CXX_FLAGS           "/GL"              )
      add_flag ( CMAKE_MODULE_LINKER_FLAGS "/LTCG:PGOPTIMIZE" )
      add_flag ( CMAKE_SHARED_LINKER_FLAGS "/LTCG:PGOPTIMIZE" )
      add_flag ( CMAKE_EXE_LINKER_FLAGS    "/LTCG:PGOPTIMIZE" )
    endif()
endif()

# common flags for intel and g++
if ( NOT WARNING_DISABLE AND ( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL ) )
   add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow" )
endif()

# C++ language features for NEC compiler
if( WALBERLA_CXX_COMPILER_IS_NEC )
   set( CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION "-Kcpp${CMAKE_CXX_STANDARD}" )
   set( CMAKE_TRY_COMPILE_PLATFORM_VARIABLES CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION )
   add_flag ( CMAKE_CXX_FLAGS "${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} -Krtti -Kexceptions -size_t64 -Kgcc" )
   add_flag ( CMAKE_CXX_FLAGS "-D__BIG_ENDIAN -D__BYTE_ORDER=__BIG_ENDIAN" )
   add_flag ( CMAKE_CXX_FLAGS "-Tnoauto,used" )
   add_flag ( CMAKE_EXE_LINKER_FLAGS "-Wl,-h,muldefs" )
   add_flag ( CMAKE_C_FLAGS "-size_t64 -Kgcc" )
   add_flag ( CMAKE_C_FLAGS "-D__BIG_ENDIAN -D__BYTE_ORDER=__BIG_ENDIAN" )
   add_flag ( CMAKE_C_FLAGS "-DSQLITE_OMIT_WAL -DHAVE_UTIME -DTHREADSAFE=0" )
   set( CMAKE_RANLIB /bin/true )
   set( CMAKE_SKIP_BUILD_RPATH TRUE )
   set( CMAKE_C_FLAGS_DEBUGOPTIMIZED    "-Chopt -g"                               )
   set( CMAKE_C_FLAGS_DEBUG             "-Cdebug -g"                              )
   set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED  "-Chopt -g"                               )
   set( CMAKE_CXX_FLAGS_DEBUG           "-Cdebug -g"                              )
endif()

# Fixes linker errors with IBM compiler
if( WALBERLA_CXX_COMPILER_IS_IBM )
   add_flag ( CMAKE_CXX_FLAGS "-qpic=large" )
endif()
# Fixes linker errors with Cray compiler
if( WALBERLA_CXX_COMPILER_IS_CRAY )
   add_flag ( CMAKE_EXE_LINKER_FLAGS  "-dynamic -L/opt/gcc/4.9.3/snos/lib64" )
endif()

# Silences compiler and linker warnings and information with the IBM compiler
if( WALBERLA_CXX_COMPILER_IS_IBM )
   add_flag ( CMAKE_CXX_FLAGS "-qsuppress=1586-267" )  # 1586-267 (I) Inlining of specified subprogram failed due to the presence of a C++ exception handler
   add_flag ( CMAKE_CXX_FLAGS "-qsuppress=1586-266" )  # 1586-266 (I) Inlining of specified subprogram failed due to the presence of a global label
   add_flag ( CMAKE_CXX_FLAGS "-qsuppress=1500-030" )  # 1500-030: (I) INFORMATION: [...] Additional optimization may be attained by recompiling and specifying MAXMEM option with a value greater than 8192.
   add_flag ( CMAKE_C_FLAGS "-qsuppress=1500-030" )    # 1500-030: (I) INFORMATION: [...] Additional optimization may be attained by recompiling and specifying MAXMEM option with a value greater than 8192.
endif()

# Silences compiler and linker warnings and information with the Cray compiler
if( WALBERLA_CXX_COMPILER_IS_CRAY )
   set( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem " )
   add_flag ( CMAKE_CXX_FLAGS "-h nomessage=1" )      # CC-1    The source file does not end with a new-line character.
   add_flag ( CMAKE_C_FLAGS   "-DSQLITE_HAVE_ISNAN" ) # SQLite will not work correctly with the -ffast-math option of GCC.
   add_flag ( CMAKE_CXX_FLAGS "-DSQLITE_HAVE_ISNAN" ) # SQLite will not work correctly with the -ffast-math option of GCC.
endif()

# architecture optimization
if( WALBERLA_OPTIMIZE_FOR_LOCALHOST )
   if( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG )
      add_flag ( CMAKE_CXX_FLAGS "-march=native" )
      add_flag ( CMAKE_C_FLAGS   "-march=native" )
      if( WALBERLA_CXX_COMPILER_IS_INTEL )
        add_flag ( CMAKE_CXX_FLAGS "-xhost" )
        add_flag ( CMAKE_C_FLAGS   "-xhost" )
      endif()
   endif()
endif()

# warning flags
if( WALBERLA_CXX_COMPILER_IS_INTEL )
   # system headers are also supported by intel, but cmake does not recognize that
   set( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem " )
   add_flag ( CMAKE_CXX_FLAGS "-wd2928,2504,2259,1682,597" )
elseif( WALBERLA_CXX_COMPILER_IS_GNU )
   add_flag ( CMAKE_CXX_FLAGS "-Wfloat-equal -Wextra" )
elseif( WALBERLA_CXX_COMPILER_IS_NEC )
   add_flag ( CMAKE_CXX_FLAGS "-wall" )
endif()

if ( WARNING_PEDANTIC AND WALBERLA_CXX_COMPILER_IS_GNU )
   add_flag ( CMAKE_CXX_FLAGS "-pedantic" )
endif ( )

 # omit deprecated warnings
if( NOT WARNING_DEPRECATED)
   if( WALBERLA_CXX_COMPILER_IS_INTEL )
       add_flag( CMAKE_CXX_FLAGS "-wd1478" )  # Disable compiler warning # 1478: "declared as deprecated"
   elseif( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_CLANG )
       add_flag ( CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")
   endif()
endif()


if ( WALBERLA_CXX_COMPILER_IS_CLANG )
    add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow -Wno-c++11-extensions -Qunused-arguments" )
endif ( )

if( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG )
    if ( WALBERLA_STL_BOUNDS_CHECKS )
        add_definitions ( "-D_GLIBCXX_DEBUG" )
        add_definitions ( "-D_LIBCPP_DEBUG=1" )
    endif()
endif()


#fastmath
if ( WALBERLA_BUILD_WITH_FASTMATH )
    if ( WALBERLA_CXX_COMPILER_IS_INTEL )
        add_flag( CMAKE_CXX_FLAGS "-fp-model fast=2 -no-prec-sqrt -no-prec-div" )
    endif()

    if ( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_CLANG )
        add_flag( CMAKE_CXX_FLAGS "-ffast-math")
    endif()

    if( WALBERLA_CXX_COMPILER_IS_MSVC )
        add_flag( CMAKE_CXX_FLAGS "/fp:fast" )
    endif()
endif()

# Xcode generator disables -isystem flag, even though current versions of Xcode support it
if(CMAKE_GENERATOR STREQUAL "Xcode")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif()


#GCC 5+ ABI selection
if( WALBERLA_CXX_COMPILER_IS_GNU )
   if( NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0.0 )
      option ( WALBERLA_USE_CPP11_ABI "On GCC 5+ use the C++11 ABI" ON )
      if( WALBERLA_USE_CPP11_ABI )
         add_flag( CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=1" )
      else()
         add_flag( CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0" )
      endif()
   endif()
endif()


# disable Xcode 7.3+ linker deduplication pass to speed up linking in debug mode
if ( APPLE )
   execute_process( COMMAND ${CMAKE_LINKER} -v OUTPUT_VARIABLE LINKER_VERSION ERROR_VARIABLE LINKER_VERSION )
   string( REGEX MATCH "ld64-[0-9\\.\\-]+" LINKER_VERSION ${LINKER_VERSION} )
   string( REGEX MATCHALL "[^\\-]+" LINKER_VERSION ${LINKER_VERSION} )
   list( GET LINKER_VERSION 0 LINKER_TYPE )
   list( GET LINKER_VERSION 1 LINKER_VERSION )
   if( LINKER_TYPE STREQUAL "ld64" AND LINKER_VERSION VERSION_GREATER 264.3.101 )
       add_flag( CMAKE_EXE_LINKER_FLAGS_DEBUG    "-Wl,-no_deduplicate")
       add_flag( CMAKE_MODULE_LINKER_FLAGS_DEBUG "-Wl,-no_deduplicate")
       add_flag( CMAKE_SHARED_LINKER_FLAGS_DEBUG "-Wl,-no_deduplicate")
   endif()
endif()


############################################################################################################################



############################################################################################################################
##
##  Find newer C++ libraries, which may only be available in std::experimental on some compilers
##
############################################################################################################################

try_compile( WALBERLA_USE_STD_FILESYSTEM "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestStdFilesystem.cpp"
             COMPILE_DEFINITIONS -DWALBERLA_USE_STD_FILESYSTEM )
if( WALBERLA_USE_STD_FILESYSTEM )
   message( STATUS "Found std::filesystem")
else()
   try_compile( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestStdFilesystem.cpp"
                COMPILE_DEFINITIONS -DWALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
   if( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
      message( STATUS "Found std::experimental::filesystem")
   endif()
   if( NOT WALBERLA_CXX_COMPILER_IS_MSVC AND NOT WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
      unset( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM CACHE )
      try_compile( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestStdFilesystem.cpp"
                   COMPILE_DEFINITIONS -DWALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM LINK_LIBRARIES stdc++fs )
      if( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
         message( STATUS "Found std::experimental::filesystem in libstdc++fs")
      list ( APPEND SERVICE_LIBS -lstdc++fs )
      endif()
   endif()
   if( NOT WALBERLA_CXX_COMPILER_IS_MSVC AND NOT WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
      unset( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM CACHE )
      try_compile( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestStdFilesystem.cpp"
                   COMPILE_DEFINITIONS -DWALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM LINK_LIBRARIES c++experimental )
      if( WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM )
         message( STATUS "Found std::experimental::filesystem in libc++experimental")
         list ( APPEND SERVICE_LIBS -lc++experimental )
      endif()
   endif()
   if( NOT WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM AND NOT WALBERLA_USE_STD_FILESYSTEM)
       message(FATAL_ERROR "Neither std::filesystem nor std::experimental::filesystem are available")
   endif()
endif()

if (WALBERLA_BUILD_WITH_CUDA AND CUDA_VERSION VERSION_LESS "11.0")
   # CUDA < 11 does not support C++17. std::experimental::any works with C++14, unlike std::any.
   set(WALBERLA_USE_STD_EXPERIMENTAL_ANY 1)
endif()


############################################################################################################################
##
##  Visual Studio Setup
##
############################################################################################################################
if ( WALBERLA_CXX_COMPILER_IS_MSVC )
   string( REGEX REPLACE "[/-]W[0-4]" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ) # remove default warning flags

   option ( WALBERLA_GROUP_PROJECTS   "Flag if the projects are grouped or in a flat hierarchy"    ON )
   option ( WALBERLA_GROUP_FILES      "Flag if the files are grouped or in a flat hierarchy"       ON )
   set_property ( GLOBAL PROPERTY USE_FOLDERS ${WALBERLA_GROUP_PROJECTS} )

   option ( WALBERLA_VS_MULTI_PROCESS_BUILD "Use the /mp option for VS builds" ON )
   if( WALBERLA_VS_MULTI_PROCESS_BUILD )
      add_flag ( CMAKE_CXX_FLAGS "-MP" ) # enable multi-threaded compiling
   endif()

   add_definitions ( "-DNOMINMAX" )                # Disable Min/Max-Macros
   add_definitions ( "-D_WIN32_WINNT=0x501" )      # Minimum Windows versions is Windows XP
   add_definitions ( "-DWINVER=0x501" )            # Minimum Windows versions is Windows XP
   add_definitions ( "-D_CRT_SECURE_NO_WARNINGS" ) # disable warnings promoting Microsoft's security enhanced CRT
   add_definitions ( "-D_SCL_SECURE_NO_WARNINGS" ) # disable warnings triggered by Microsoft's checked iterators
   add_flag ( CMAKE_CXX_FLAGS "-W4" )              # set warning level to maximum
   add_flag ( CMAKE_CXX_FLAGS "-bigobj" )          # enable big object files
   add_flag ( CMAKE_CXX_FLAGS "-wd4127" )          # disable compiler warning C4127: "conditional expression is constant"
   add_flag ( CMAKE_CXX_FLAGS "-wd4512" )          # disable compiler warning C4512: "assignment operator could not be generated"
   add_flag ( CMAKE_CXX_FLAGS "-wd4913" )          # disable compiler warning C4512: "user defined binary operator ',' exists but
                                                   # no overload could convert all operands, default built-in binary operator ','
                                                   # used"
   add_flag ( CMAKE_CXX_FLAGS "-wd4702" )          # disable compiler warning C4702: "unreachable code"
   add_flag ( CMAKE_CXX_FLAGS "-wd4505" )          # disable compiler warning C4505: "unreferenced local function has been removed"
   add_flag ( CMAKE_CXX_FLAGS "-wd4503" )          # disable compiler warning C4503: "'identifier' : decorated name length exceeded, name was truncated"

   if ( WARNING_ERROR )
      add_flag ( CMAKE_CXX_FLAGS "-WX" )           # Treat warnings as errors
   endif ( )

   if( NOT WARNING_DEPRECATED)
      add_definitions( "-D_CRT_SECURE_NO_DEPRECATE" )
      add_definitions( "-D_SCL_SECURE_NO_DEPRECATE" )
      add_flag       ( CMAKE_CXX_FLAGS "-wd4996"    ) # Disable compiler warning C4996: "declared as deprecated"
   endif()

endif ( )
############################################################################################################################




############################################################################################################################
##
## Code Generation (pystencils)
##
#############################################################################################################################
if ( WALBERLA_BUILD_WITH_CODEGEN )
   find_package( PythonInterp 3 QUIET REQUIRED)
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import lbmpy"
         RESULT_VARIABLE LBMPY_FOUND )
    if(NOT LBMPY_FOUND EQUAL 0)
       message(FATAL_ERROR "WALBERLA_BUILD_WITH_CODEGEN activated but pystencils or lbmpy package not found.
                            Please install lbmpy e.g.: 'pip3 install lbmpy'")
    endif()
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from pystencils.include import get_pystencils_include_path; print(get_pystencils_include_path())"
                    OUTPUT_VARIABLE PYSTENCILS_INCLUDE_PATH)
    include_directories( ${PYSTENCILS_INCLUDE_PATH} )

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import jinja2"
          RESULT_VARIABLE JINJA2_FOUND )
    if(NOT JINJA2_FOUND EQUAL 0)
       message(FATAL_ERROR "WALBERLA_BUILD_WITH_CODEGEN activated and jinja2 package not found.
       Please install jinja2 e.g.: 'pip3 install jinja2'")
    endif()
endif()
############################################################################################################################




############################################################################################################################
##
## Python Libraries
##
#############################################################################################################################
if ( WALBERLA_BUILD_WITH_PYTHON )

    find_package( PythonInterp 3.6 QUIET REQUIRED)
    find_package( PythonLibs QUIET REQUIRED)

    if( NOT (PYTHON_LIBRARY AND PYTHON_INCLUDE_DIR ) )
        message( FATAL_ERROR "Couldn't find any python library" )
    endif()

    if(WALBERLA_CXX_COMPILER_IS_INTEL)
        # Intel C++17 support introduced in 2.6.2 (https://github.com/pybind/pybind11/pull/2729)
        set(PYBIND11_MINIMUM_VERSION "2.6.2")
    else()
        set(PYBIND11_MINIMUM_VERSION "2.6.0")
    endif()

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pybind11; print(pybind11._version.__version__)"
                    OUTPUT_VARIABLE pybind11_VERSION ERROR_QUIET RESULT_VARIABLE pybind11_VERSION_RESULT)
    string(STRIP "${pybind11_VERSION}" pybind11_VERSION)
    if(pybind11_VERSION_RESULT EQUAL "0" AND pybind11_VERSION VERSION_GREATER_EQUAL "${PYBIND11_MINIMUM_VERSION}")
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
                        OUTPUT_VARIABLE PYBIND11_CMAKE_PATH)
        string(STRIP "${PYBIND11_CMAKE_PATH}" PYBIND11_CMAKE_PATH)
        find_package(pybind11 PATHS "${PYBIND11_CMAKE_PATH}" NO_DEFAULT_PATH REQUIRED)
    else()
        find_package(Git QUIET)
        if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
           # Update submodules as needed
           if(WALBERLA_GIT_SUBMODULE_AUTO)
              message(STATUS "Submodule update")
              execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
              if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                 message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
              endif()
           endif()
        endif()

        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt")
            if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
                message(FATAL_ERROR "Please update git submodules or install pybind11 via pip.")
            else()
                message(FATAL_ERROR "Please install pybind11 via pip or download it to ${PROJECT_SOURCE_DIR}/extern/pybind11.")
            endif()
        endif()

        add_subdirectory(extern/pybind11)

        if(pybind11_VERSION VERSION_LESS "2.6.2")
            # if pybind11 was installed into ${PYTHON_INCLUDE_DIR} (e.g. by pip), that will have a higher priority than the Git submodule unless we reorder the search path
            # introduced in 2.6.0 (https://github.com/pybind/pybind11/issues/2709), fixed in 2.6.2 (https://github.com/pybind/pybind11/pull/2716)
            set_property( TARGET pybind11::pybind11
                          PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/extern/pybind11/include" "${PYTHON_INCLUDE_DIR}")
        endif()
    endif()

    # a python module is a shared library - so everything has to be compiled to position independent code
    # otherwise linking the static libs into the shared lib will result in errors
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)


    if(WALBERLA_BUILD_DOC)
      # Sphinx documentation
      # to build documentation make sure to have sphinx and read-the-docs theme installed
      # Install with: "pip install sphinx sphinx_rtd_theme"
      add_custom_target( docPython sphinx-build -b html "${walberla_SOURCE_DIR}/python/waLBerla_docs" "${walberla_BINARY_DIR}/doc/python"
                         COMMENT "Building HTML documentation for Python extension with Sphinx")
  endif()
endif()


############################################################################################################################
##
## BOOST Libraries
##
#############################################################################################################################
set ( Boost_NO_BOOST_CMAKE ON)

# This variable is necessary, if the CMAKE version used is not aware of a more recent boost version (keep this up to date!)
set ( Boost_ADDITIONAL_VERSIONS
      "1.45" "1.45.0" "1.46" "1.46.0" "1.46.1" "1.47" "1.47.0" "1.48" "1.48.0" "1.49" "1.49.0"
      "1.50" "1.50.0" "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0" "1.54" "1.54.0" "1.55" "1.55.0"
      "1.56" "1.56.0" "1.57" "1.57.0" "1.58" "1.58.0" "1.59" "1.59.0" "1.60" "1.60.0" "1.61" "1.61.0" "1.62" "1.62.0" "1.63" "1.63.0"
      "1.64.0" "1.65.0" "1.65.1" "1.66.0" "1.67.0" "1.68.0" "1.69.0" "1.70.0" "1.71.0" "1.72.0" "1.73.0" "1.74.0")

# if you defined BOOST_ROOT or BOOST_BASE in your environment use it here to find boost too
if ( NOT BOOST_ROOT )
   foreach ( var  BOOST_ROOT  BOOST_BASE )
      if ( NOT "$ENV{${var}}" STREQUAL "" )
         message ( STATUS "Use environment boost directory: $ENV{${var}}" )
         set ( BOOST_ROOT $ENV{${var}} CACHE INTERNAL "")
         break ( )
      endif ( )
   endforeach ( )
endif ( )

find_package ( Boost )

if ( Boost_FOUND )
   if(CMAKE_GENERATOR STREQUAL "Xcode")
      # this is needed because the SYSTEM flag to include_directories does not work
      add_flag ( CMAKE_CXX_FLAGS "-isystem ${Boost_INCLUDE_DIRS}" )
   else()
      include_directories ( SYSTEM ${Boost_INCLUDE_DIRS} )
   endif()
   add_definitions ( -DBOOST_ALL_NO_LIB ) # Disable Boost auto-linking (CMAKE does that for us...)

   set( WALBERLA_BUILD_WITH_BOOST TRUE CACHE INTERNAL "Build with Boost" )
else( Boost_FOUND )
   set( WALBERLA_BUILD_WITH_BOOST FALSE CACHE INTERNAL "Build with Boost" )
endif( Boost_FOUND )


############################################################################################################################





############################################################################################################################
##
## PThread is required in Linux environments by std::thread
##
############################################################################################################################

if ( NOT WIN32 )
   add_flag( CMAKE_CXX_FLAGS "-pthread" )
endif()



############################################################################################################################
##
## backtrace may be in a separate library
##
############################################################################################################################

if ( NOT WIN32 AND (WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG))
   find_package ( Backtrace QUIET )
   if ( Backtrace_FOUND )
      list ( APPEND SERVICE_LIBS ${Backtrace_LIBRARIES} )
      set ( WALBERLA_BUILD_WITH_BACKTRACE ON )
      set ( WALBERLA_BACKTRACE_HEADER ${Backtrace_HEADER} )
   endif ( Backtrace_FOUND )
endif()



############################################################################################################################
##
## MPI
##
############################################################################################################################

if ( WALBERLA_BUILD_WITH_MPI AND NOT WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER )
   find_package ( MPI )

   # FindMPI does not really work under windows, because it expects linux formatted strings from the mpi compiler.
   # Nevertheless for Microsoft MPI and MPICH there are workarounds included, but not for OpenMPI.
   # Here is a workaround for windows with OpenMPI (use configure twice to see correct output).
   # The workaround bases on the elseif(try_libs) case in the interrogate_* function of FindMPI.
   # For this workaround we assume, that the compiler was found in any PATH or ENV variable
   # but FindMPI was not able to interpret the command line outputs.
   if ( MPI_CXX_COMPILER AND NOT MPI_CXX_LIBRARIES )
      if ( WIN32 )
         message ( STATUS "Enter Workaround Routine for Windows and OpenMPI: PRESS CONFIGURE ONE MORE TIME!" )
         string ( REGEX REPLACE "(.*)/bin/.*" "\\1" MPI_PATH ${MPI_CXX_COMPILER} )
         find_path ( MPI_C_INCLUDE_PATH mpi.h
            HINTS ${MPI_PATH}
            PATH_SUFFIXES include Inc)
         set ( MPI_CXX_INCLUDE_PATH ${MPI_C_INCLUDE_PATH} CACHE FILEPATH "" FORCE )

         set ( MPI_CXX_LIBRARIES "MPI_CXX_LIBRARIES-NOTFOUND" CACHE FILEPATH "Cleared" FORCE )
         find_library ( MPI_CXX_LIBRARIES
            NAMES         mpi++ mpicxx cxx mpi_cxx libmpi++ libmpicxx libcxx libmpi_cxx
            HINTS         ${MPI_PATH}
            PATH_SUFFIXES lib )

         if ( NOT MPI_CXX_LIBRARIES STREQUAL "MPI_CXX_LIBRARIES-NOTFOUND" )
            set ( MPI_CXX_FOUND ON FORCE )
         endif ( )

         set ( MPI_C_LIBRARIES "MPI_C_LIBRARIES-NOTFOUND" CACHE FILEPATH "Cleared" FORCE )
         find_library ( MPI_C_LIBRARIES
           NAMES         mpi mpich mpich2 msmpi libmpi libmpich libmpich2 libmsmpi
           HINTS         ${MPI_PATH}
           PATH_SUFFIXES lib )

         if ( NOT MPI_C_LIBRARIES STREQUAL "MPI_C_LIBRARIES-NOTFOUND" )
            set ( MPI_C_FOUND ON FORCE )
         endif ( )

         if ( MPI_PATH MATCHES ".*OpenMPI.*" )
            set ( MPI_CXX_COMPILE_FLAGS "/DOMPI_IMPORTS" CACHE STRING "" FORCE )
            set ( MPI_C_COMPILE_FLAGS   "/DOMPI_IMPORTS" CACHE STRING "" FORCE )
         endif ( )
      elseif ( WALBERLA_CXX_COMPILER_IS_CRAY )
      else ( )
         message ( WARNING "Found MPI Compiler but no Libraries -> invent a new workaround" )
      endif ( )
   endif ( )

   if ( MPI_FOUND )
     include_directories ( SYSTEM ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH} )
     foreach( LIB ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} )
         if ( LIB )
            list ( APPEND SERVICE_LIBS ${LIB} )
         endif ( )
     endforeach ( )
     add_flag ( CMAKE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}" )
     add_flag ( CMAKE_C_FLAGS   "${MPI_C_COMPILE_FLAGS}" )

     add_flag ( CMAKE_MODULE_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}" )
     add_flag ( CMAKE_EXE_LINKER_FLAGS    "${MPI_CXX_LINK_FLAGS}" )
     add_flag ( CMAKE_SHARED_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}" )
   endif ( )
endif ( )

# OpenMPI 3.0 and higher checks the number of processes against the number of CPUs
execute_process(COMMAND ${MPIEXEC} --version RESULT_VARIABLE mpi_version_result OUTPUT_VARIABLE mpi_version_output)
if (mpi_version_result EQUAL 0 AND mpi_version_output MATCHES "\\(Open(RTE| MPI)\\) ([3-9]\\.|1[0-9])")
    set ( MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS}" "-oversubscribe" CACHE STRING "" FORCE)
endif()
############################################################################################################################



############################################################################################################################
##
## Qt
##
############################################################################################################################
option (WALBERLA_ENABLE_GUI "This flag builds the graphical user interface, depends on Qt Libraries")

if ( WALBERLA_ENABLE_GUI )

    find_package( Qt4 COMPONENTS QtCore QtGui QtOpenGL QtXml REQUIRED )
    set( OpenGL_GL_PREFERENCE LEGACY )
    find_package( OpenGL REQUIRED )

    INCLUDE( ${QT_USE_FILE} )
    list ( APPEND SERVICE_LIBS ${OPENGL_LIBRARIES} ${QT_LIBRARIES} )

endif(WALBERLA_ENABLE_GUI)

############################################################################################################################



############################################################################################################################
##
## METIS
##
############################################################################################################################

if ( WALBERLA_BUILD_WITH_METIS )
    find_package ( Metis QUIET )

    if ( METIS_FOUND )
        include_directories( ${METIS_INCLUDE_DIRS} )
        link_directories   ( ${METIS_LIBRARY_DIR}  )
        list ( APPEND SERVICE_LIBS ${METIS_LIBRARIES} )
        set  ( WALBERLA_BUILD_WITH_METIS TRUE )
    else()
        set  ( WALBERLA_BUILD_WITH_METIS OFF CACHE BOOL "Build with metis graph partitioner" FORCE )
    endif()
else()
    set ( METIS_FOUND OFF CACHE BOOL "Metis found" FORCE )
endif()


if ( WALBERLA_BUILD_WITH_PARMETIS )
   find_path(PARMETIS_INCLUDE_DIR parmetis.h
      /usr/local/include
      /usr/include
      ${PARMETIS_ROOT}/include
      $ENV{PARMETIS_ROOT}/include
   )

  find_library(PARMETIS_LIBRARY parmetis
    /usr/local/lib
    /usr/lib
    ${PARMETIS_ROOT}/lib
    $ENV{PARMETIS_ROOT}/lib
  )

  if( PARMETIS_INCLUDE_DIR AND PARMETIS_LIBRARY AND METIS_LIBRARY )
    include_directories( ${PARMETIS_INCLUDE_DIR} )
    list ( APPEND SERVICE_LIBS ${PARMETIS_LIBRARY} ${METIS_LIBRARY} )
  endif()
endif()

############################################################################################################################



############################################################################################################################
##
## FFTW3
##
############################################################################################################################

if( WALBERLA_BUILD_WITH_MPI )
   find_package( PFFT )
   find_package( FFTW3 )
   set( FFT_REQUIRED_LIBRARIES pfft fftw3_mpi fftw3 )
   if( PFFT_FOUND AND FFTW3_MPI_FOUND )
      set( WALBERLA_BUILD_WITH_FFT TRUE CACHE INTERNAL "Build with FFT" )
      include_directories( SYSTEM ${PFFT_INCLUDE_DIR} ${FFTW3_MPI_INCLUDE_DIR} )
      list( APPEND SERVICE_LIBS ${PFFT_LIBRARIES} ${FFTW3_LIBRARIES} ${FFTW3_MPI_LIBRARIES} )
   endif()
else()
   find_package( FFTW3 )
   set( FFT_REQUIRED_LIBRARIES fftw3 )
   if ( FFTW3_FOUND )
      set( WALBERLA_BUILD_WITH_FFT TRUE CACHE INTERNAL "Build with FFT" )
      include_directories( SYSTEM ${FFTW3_INCLUDE_DIR} )
      list( APPEND SERVICE_LIBS ${FFTW3_LIBRARIES} )
   endif()
endif()



############################################################################################################################
##
## OpenMesh
##
############################################################################################################################
if( (NOT DEFINED WALBERLA_BUILD_WITH_OPENMESH) OR WALBERLA_BUILD_WITH_OPENMESH )
   find_package( OpenMesh )
   if( OPENMESH_FOUND )
      set( WALBERLA_BUILD_WITH_OPENMESH ON CACHE BOOL "Build with OpenMesh support" )
      include_directories( SYSTEM ${OPENMESH_INCLUDE_DIRS} )
      list( APPEND SERVICE_LIBS ${OPENMESH_LIBRARIES} )
      if( WALBERLA_CXX_COMPILER_IS_MSVC )
         add_definitions(-D_USE_MATH_DEFINES )
      endif()
   else()
      set( WALBERLA_BUILD_WITH_OPENMESH OFF CACHE BOOL "Build with OpenMesh support" FORCE )
   endif()
endif()

############################################################################################################################
##
## DebugOptimized Build Configuration for fast execution of tests with enabled asserts
##
############################################################################################################################

set( CMAKE_C_FLAGS_DEBUGOPTIMIZED             ${CMAKE_C_FLAGS_DEBUG}             )
set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED           ${CMAKE_CXX_FLAGS_DEBUG}           )
set( CMAKE_EXE_LINKER_FLAGS_DEBUGOPTIMIZED    ${CMAKE_EXE_LINKER_FLAGS_DEBUG}    )
set( CMAKE_SHARED_LINKER_FLAGS_DEBUGOPTIMIZED ${CMAKE_SHARED_LINKER_FLAGS_DEBUG} )
set( CMAKE_MODULE_LINKER_FLAGS_DEBUGOPTIMIZED ${CMAKE_MODULE_LINKER_FLAGS_DEBUG} )

set_property(GLOBAL PROPERTY DEBUG_CONFIGURATIONS Debug DebugOptimized)

if ( WALBERLA_CXX_COMPILER_IS_MSVC )
    string(REPLACE "/Od" "/O2"   CMAKE_C_FLAGS_DEBUGOPTIMIZED   ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
    string(REPLACE "/Ob0" "/Ob2" CMAKE_C_FLAGS_DEBUGOPTIMIZED   ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
    string(REPLACE "/RTC1" ""    CMAKE_C_FLAGS_DEBUGOPTIMIZED   ${CMAKE_C_FLAGS_DEBUGOPTIMIZED})
    string(REPLACE "/Od" "/O2"   CMAKE_CXX_FLAGS_DEBUGOPTIMIZED ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
    string(REPLACE "/Ob0" "/Ob2" CMAKE_CXX_FLAGS_DEBUGOPTIMIZED ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
    string(REPLACE "/RTC1" ""    CMAKE_CXX_FLAGS_DEBUGOPTIMIZED ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED})
elseif( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG )
   set( CMAKE_C_FLAGS_DEBUGOPTIMIZED   "${CMAKE_C_FLAGS_DEBUGOPTIMIZED} -O3" )
   set( CMAKE_CXX_FLAGS_DEBUGOPTIMIZED "${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED} -O3" )
endif()

set(CMAKE_C_FLAGS_DEBUGOPTIMIZED ${CMAKE_C_FLAGS_DEBUGOPTIMIZED} CACHE STRING
    "Flags used by the compiler during DebugOptimized builds")
set(CMAKE_CXX_FLAGS_DEBUGOPTIMIZED ${CMAKE_CXX_FLAGS_DEBUGOPTIMIZED}  CACHE STRING
    "Flags used by the compiler during DebugOptimized builds")
set(CMAKE_EXE_LINKER_FLAGS_DEBUGOPTIMIZED ${CMAKE_EXE_LINKER_FLAGS_DEBUGOPTIMIZED} CACHE STRING
    "Flags used by the linker for executables during DebugOptimized builds")
set(CMAKE_SHARED_LINKER_FLAGS_DEBUGOPTIMIZED ${CMAKE_SHARED_LINKER_FLAGS_DEBUGOPTIMIZED} CACHE STRING
    "Flags used by the linker for shared libraries during DebugOptimized builds")
set(CMAKE_MODULE_LINKER_FLAGS_DEBUGOPTIMIZED ${CMAKE_MODULE_LINKER_FLAGS_DEBUGOPTIMIZED} CACHE STRING
    "Flags used by the linker for loadable modules during DebugOptimized builds")

mark_as_advanced(
    CMAKE_C_FLAGS_DEBUGOPTIMIZED
    CMAKE_CXX_FLAGS_DEBUGOPTIMIZED
    CMAKE_EXE_LINKER_FLAGS_DEBUGOPTIMIZED
    CMAKE_SHARED_LINKER_FLAGS_DEBUGOPTIMIZED
    CMAKE_MODULE_LINKER_FLAGS_DEBUGOPTIMIZED
)

# This variable is only set for multi-config IDE generators like VC
if( CMAKE_CONFIGURATION_TYPES )
    list( APPEND CMAKE_CONFIGURATION_TYPES DebugOptimized )
    list( REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES )
    set ( CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}"
          CACHE STRING "Semicolon separated list of supported types
                      [Debug|Release|MinSizeRel|RelWithDebInfo|ReleaseNoOutfiles|DebugOptimized]" FORCE )
endif()
############################################################################################################################





############################################################################################################################
##
## OpenMP
##
############################################################################################################################

option ( WALBERLA_THREAD_SAFE_LOGGING "Enables/Disables thread-safe logging" ON )

if ( WALBERLA_BUILD_WITH_OPENMP )
    if ( WALBERLA_CXX_COMPILER_IS_INTEL AND "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "16.0.3" )
       add_flag ( CMAKE_C_FLAGS   "-openmp" )
       add_flag ( CMAKE_CXX_FLAGS "-openmp" )
    elseif ( WALBERLA_CXX_COMPILER_IS_NEC )
       add_flag ( CMAKE_C_FLAGS   "-Popenmp" )
       add_flag ( CMAKE_CXX_FLAGS "-Popenmp" )
    else()
       find_package( OpenMP )
       add_flag ( CMAKE_C_FLAGS   "${OpenMP_C_FLAGS}" )
       add_flag ( CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}" )
       list ( APPEND SERVICE_LIBS ${OpenMP_CXX_LIBRARIES} )
    endif()
else()
    if ( WALBERLA_CXX_COMPILER_IS_CRAY )
       add_flag ( CMAKE_C_FLAGS   "-h noomp" )
       add_flag ( CMAKE_CXX_FLAGS "-h noomp" )
    endif()
endif()
############################################################################################################################



############################################################################################################################
##
## CUDA
##
############################################################################################################################
if ( WALBERLA_BUILD_WITH_CUDA )
    include(CheckLanguage)
    check_language(CUDA)
    if( CMAKE_CUDA_COMPILER )
      enable_language(CUDA)

      #include directories and cudart lib is needed for cpp files that use cuda headers/libs
      include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
      find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
      list ( APPEND SERVICE_LIBS ${CUDART_LIBRARY} )

      #TODO: check if this is really needed? is it possible that nvtx is missing?
      find_library( NVTX_LIBRARY nvToolsExt ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES} )
      if( NVTX_LIBRARY )
          set( WALBERLA_CUDA_NVTX_AVAILABLE 1)
          list ( APPEND SERVICE_LIBS ${NVTX_LIBRARY} )
      endif()

      #CUDA_FOUND is need for our cmake mechanism
      set ( CUDA_FOUND TRUE )
    else()
      set ( WALBERLA_BUILD_WITH_CUDA FALSE )
    endif ( )
endif ( )
############################################################################################################################



############################################################################################################################
##
##  Testing Coverage
##
############################################################################################################################
if (WALBERLA_BUILD_WITH_GCOV AND CMAKE_COMPILER_IS_GNUCXX  )
    add_flag ( CMAKE_CXX_FLAGS "--coverage" )
endif()
############################################################################################################################



############################################################################################################################
##
##  Profiling with gprof
##
############################################################################################################################

if ( WALBERLA_BUILD_WITH_GPROF )
    if ( WALBERLA_CXX_COMPILER_IS_INTEL )
        add_flag ( CMAKE_CXX_FLAGS        "-pg" )
        add_flag ( CMAKE_EXE_LINKER_FLAGS "-pg" )
    elseif ( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_CLANG )
        add_flag ( CMAKE_CXX_FLAGS        "-pg" )
    endif()
endif()
############################################################################################################################


############################################################################################################################
##
##  Likwid Marker API
##
############################################################################################################################


if ( WALBERLA_BUILD_WITH_LIKWID_MARKERS )
    find_library( LIKWID_LIB likwid HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )
    find_path( LIKWID_INCLUDE_DIR likwid.h HINTS $ENV{LIKWID_INCDIR} $ENV{LIKWID_ROOT}/include )

    if ( LIKWID_LIB AND LIKWID_INCLUDE_DIR)
        set( LIKWID_FOUND 1)
        include_directories( ${LIKWID_INCLUDE_DIR})
        add_definitions ( "-DLIKWID_PERFMON" )
        list ( APPEND SERVICE_LIBS ${LIKWID_LIB} )
    else()
        message(WARNING "likwid marker library not found. Set environment variable LIKWID_ROOT")
        set ( WALBERLA_BUILD_WITH_LIKWID_MARKERS OFF CACHE BOOL "Compile in markers for likwid-perfctr" FORCE )
    endif()
endif()



############################################################################################################################
##
##  Link time optimizations
##
############################################################################################################################
if ( WALBERLA_BUILD_WITH_LTO  )

   if( WALBERLA_CXX_COMPILER_IS_INTEL )
      add_flag( CMAKE_CXX_FLAGS_RELEASE "-ip -ipo3" )
      add_flag( CMAKE_C_FLAGS_RELEASE   "-ip -ipo3" )
   endif()

   if ( CMAKE_COMPILER_IS_GNUCXX )
      add_flag ( CMAKE_C_FLAGS_RELEASE     "-flto=3" )
      add_flag ( CMAKE_CXX_FLAGS_RELEASE   "-flto=3" )
      add_flag ( CMAKE_EXE_LINKER_FLAGS    "-fuse-linker-plugin" )
   endif ( )

   if( WALBERLA_CXX_COMPILER_IS_MSVC )
      add_flag ( CMAKE_CXX_FLAGS_RELEASE           "/GL"   )
      add_flag ( CMAKE_EXE_LINKER_FLAGS_RELEASE    "/LTCG" )
      add_flag ( CMAKE_SHARED_LINKER_FLAGS_RELEASE "/LTCG" )
      add_flag ( CMAKE_MODULE_LINKER_FLAGS_RELEASE "/LTCG" )
   endif ( )

   if( WALBERLA_CXX_COMPILER_IS_IBM )
      add_flag ( CMAKE_C_FLAGS_RELEASE     "-qipa" )
      add_flag ( CMAKE_CXX_FLAGS_RELEASE   "-qipa" )
      add_flag ( CMAKE_EXE_LINKER_FLAGS    "-qipa" )
   endif( )

endif ( )

############################################################################################################################
##
##  Some more compiler flags that need to happen after any try_compile (e.g. inside FindMPI)
##
############################################################################################################################

# Treat warnings as errors
if ( WARNING_ERROR )
   if( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG )
      add_flag ( CMAKE_CXX_FLAGS "-pedantic-errors -Werror" )
   elseif( WALBERLA_CXX_COMPILER_IS_MSVC )
      add_flag ( CMAKE_CXX_FLAGS "/WX" )
   elseif ( WALBERLA_CXX_COMPILER_IS_CRAY )
      add_flag ( CMAKE_CXX_FLAGS "-h error_on_warning" )
   endif()
endif ( )

############################################################################################################################
##
##  Sanitizer
##
############################################################################################################################
if ( WALBERLA_SANITIZE_ADDRESS )
    if ( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_CLANG )
        add_flag( CMAKE_CXX_FLAGS "-fsanitize=address")
    endif()
endif()

if ( WALBERLA_SANITIZE_UNDEFINED )
    if ( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_CLANG )
        add_flag( CMAKE_CXX_FLAGS "-fsanitize=undefined")
    endif()
endif()

############################################################################################################################
# Documentation Generation
#
if (WALBERLA_BUILD_DOC)
  # Build documentation using Doxygen (www.doxygen.org)
  ############################################################################################################################
  find_package ( Doxygen  )
  find_package ( HTMLHelp )

  if ( HTML_HELP_COMPILER EQUAL "" )
     set ( HTML_HELP_FOUND "NO" )
  else ( )
     set ( HTML_HELP_FOUND "YES" )
  endif ( )

  if ( DOXYGEN_FOUND )
     set ( DOXYGEN_HTML_HEADER ${walberla_SOURCE_DIR}/doc/header.html )
     set ( DOXYGEN_HTML_FOOTER ${walberla_SOURCE_DIR}/doc/footer.html )
     set ( DOXYGEN_HTML_OUTPUT "html" )

     configure_file ( ${walberla_SOURCE_DIR}/doc/doxygen.in ${walberla_BINARY_DIR}/doc/doxygen.cfg @ONLY )

     add_custom_target ( doc   ${DOXYGEN_EXECUTABLE} ${walberla_BINARY_DIR}/doc/doxygen.cfg
                            COMMENT "Generating API documentation with Doxygen" VERBATIM )

  endif ( )
endif()

############################################################################################################################





############################################################################################################################
#
# Fix compiler bugs
#
############################################################################################################################

# The NEC SX has a few issues in its standard library headers
if( WALBERLA_CXX_COMPILER_IS_NEC )
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/math.h         "#include_next <math.h>\n#undef fpclassify\n#undef signbit\n#undef isfinite\n#undef isinf\n#undef isnan\n#undef isnormal\n#undef isgreater\n#undef isgreaterequal\n#undef isless\n#undef islessequal\n#undef islessgreater\n#undef isunordered\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/types.h    "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/types.h\"   \n#undef uint_t\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/acl.h      "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/acl.h\"     \n#undef uint_t\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/if_ehcpl.h "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/if_ehcpl.h\"\n#undef uint_t\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/ptms.h     "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/ptms.h\"    \n#undef uint_t\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/stream.h   "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/stream.h\"  \n#undef uint_t\n")
   file( WRITE ${walberla_BINARY_DIR}/CMakeFiles/src/sys/strsubr.h  "#define uint_t SX_UINT_T\n#include \"/SX/usr/include/sys/strsubr.h\" \n#undef uint_t\n")
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/math.h         ${walberla_BINARY_DIR}/src/math.h COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/types.h    ${walberla_BINARY_DIR}/src/sys/types.h    COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/acl.h      ${walberla_BINARY_DIR}/src/sys/acl.h      COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/if_ehcpl.h ${walberla_BINARY_DIR}/src/sys/if_ehcpl.h COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/ptms.h     ${walberla_BINARY_DIR}/src/sys/ptms.h     COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/stream.h   ${walberla_BINARY_DIR}/src/sys/stream.h   COPYONLY )
   configure_file ( ${walberla_BINARY_DIR}/CMakeFiles/src/sys/strsubr.h  ${walberla_BINARY_DIR}/src/sys/strsubr.h  COPYONLY )
endif()



############################################################################################################################
#
# Subdirectories
#
############################################################################################################################

# Add binary dir for generated headers
include_directories ( ${CMAKE_CURRENT_BINARY_DIR}/src )
# All include paths are specified relative to src/ directory
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/src )


# Generate file with compile options, and add install rule for it
configure_file ( src/waLBerlaDefinitions.in.h
                 src/waLBerlaDefinitions.h    )

install( FILES ${walberla_BINARY_DIR}/src/waLBerlaDefinitions.h DESTINATION walberla/ )

# sources
add_subdirectory ( src )

# test
if ( WALBERLA_BUILD_TESTS )
    add_subdirectory ( tests )
else()
    add_subdirectory( tests EXCLUDE_FROM_ALL )
endif()


add_subdirectory ( apps )

waLBerla_export()

############################################################################################################################

############################################################################################################################
##
## clang-tidy
##
############################################################################################################################

waLBerla_link_files_to_builddir( .clang-tidy )
add_subdirectory( utilities )
