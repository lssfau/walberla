############################################################################################################################
##
## waLBerla's main cmake file
##
## Contents:
##   - definition of build options
##   - compiler variables ( c++ standard, warnings etc. )
##   - Finding of service libraries. Required: none, Optional: MPI, FFTW3, METIS, OpenMesh, Python
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

CMAKE_MINIMUM_REQUIRED (VERSION 3.24)


PROJECT ( walberla )

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${walberla_SOURCE_DIR}/cmake ${walberla_SOURCE_DIR}/cmake/compileroptions )

include ( waLBerlaFunctions )

set_version( 7 1 )

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

option ( WALBERLA_BUILD_TESTS               "Build Testcases"                                 ON )
option ( WALBERLA_BUILD_BENCHMARKS          "Build Benchmarks"                                ON )
option ( WALBERLA_BUILD_TOOLS               "Build Tools"                                     ON )
option ( WALBERLA_BUILD_TUTORIALS           "Build Tutorials"                                 ON )
option ( WALBERLA_BUILD_SHOWCASES           "Build Showcases"                                 ON )
option ( WALBERLA_BUILD_DOC                 "Build Documentation"                             ON )

option ( WALBERLA_BUILD_WITH_MPI            "Build with MPI"                                  ON )
option ( WALBERLA_BUILD_WITH_METIS          "Build with metis graph partitioner"             OFF )
option ( WALBERLA_BUILD_WITH_PARMETIS       "Build with ParMetis graph partitioner"          OFF )
option ( WALBERLA_BUILD_WITH_FFTW           "Build with FFTW Fourier Transform library"      OFF )

option ( WALBERLA_BUILD_WITH_GPROF          "Enables gprof"                                      )
option ( WALBERLA_BUILD_WITH_GCOV           "Enables gcov"                                       )
option ( WALBERLA_BUILD_WITH_LTO            "Enable link time optimizations"                     )
option ( WALBERLA_BUILD_WITH_OPENMP         "Enable OpenMP support"                              )
option ( WALBERLA_BUILD_WITH_PYTHON         "Support for embedding Python"                       )
option ( WALBERLA_BUILD_WITH_CODEGEN        "Enable pystencils code generation"              OFF )


option ( WALBERLA_BUILD_WITH_LIKWID_MARKERS "Compile in markers for likwid-perfctr"              )

option ( WALBERLA_BUILD_WITH_CUDA	        "Enable CUDA support"                                )
option ( WALBERLA_BUILD_WITH_HIP	           "Enable ROCm HIP support"                            )


option ( WALBERLA_BUILD_WITH_FASTMATH       "Fast math"                                          )

option ( WALBERLA_SIMD_FORCE_SCALAR         "Do not use SIMD operations even when available" OFF )

option ( WALBERLA_BUFFER_DEBUG              "Type checking for BufferSystem ( slow )"        OFF )

# Profile guided optimization
option ( WALBERLA_PROFILE_GENERATE          "Generates Profile for Optimization"                 )
option ( WALBERLA_PROFILE_USE               "Uses Profile to optimize"                           )

# Compiler Optimization
option ( WALBERLA_OPTIMIZE_FOR_LOCALHOST    "Enable compiler optimizations spcific to localhost" )

option ( WALBERLA_LOG_SKIPPED               "Log skipped cmake targets"                      ON  )
option ( WALBERLA_DEPS_ERROR                "Fail if module dependencies are not met"        OFF )

option ( WALBERLA_GIT_SUBMODULE_AUTO        "Check submodules during cmake run"               ON )

option ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT "Experimental half precision support"    OFF )

# Installation Directory
set ( CMAKE_INSTALL_PREFIX /usr/local/waLBerla CACHE STRING "The default installation directory."   )

# Default build type
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release CACHE STRING "Build Types: Debug Release DebugOptimized RelWithDebInfo MinSizeRel."  FORCE )
endif()
SET_PROPERTY( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release DebugOptimized RelWithDebInfo MinSizeRel )
set( CMAKE_TRY_COMPILE_CONFIGURATION ${CMAKE_BUILD_TYPE} )

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
     include(GNU)
else()
     option ( WALBERLA_CXX_COMPILER_IS_GNU "Use gnu compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_GNU )

# Check for Visual Studio
if ( MSVC )
     option ( WALBERLA_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" ON  )
     include( MSVC )
else()
     option ( WALBERLA_CXX_COMPILER_IS_MSVC "Use Visual Studio compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_MSVC )

# Check for IBM compiler
if( CMAKE_CXX_COMPILER MATCHES "xlc" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "xlc" )
    option ( WALBERLA_CXX_COMPILER_IS_IBM "Use IBM compiler" ON  )
    include(IBM)
else()
    option ( WALBERLA_CXX_COMPILER_IS_IBM "Use IBM compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_IBM )

if( CMAKE_CXX_COMPILER MATCHES "clang" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "clang" OR CMAKE_CXX_COMPILER MATCHES "hipcc" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" )
    option ( WALBERLA_CXX_COMPILER_IS_CLANG "Use clang compiler" ON  )
    include(Clang)
else()
    option ( WALBERLA_CXX_COMPILER_IS_CLANG "Use clang compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_CLANG )

# Check for Cray compiler
if( CMAKE_CXX_COMPILER_ID MATCHES Cray )
    option ( WALBERLA_CXX_COMPILER_IS_CRAY "Use Cray compiler" ON   )
    include(Cray)
else()
    option ( WALBERLA_CXX_COMPILER_IS_CRAY "Use Cray compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_CRAY )

# Check for Fujitsu compiler
if( CMAKE_CXX_COMPILER_ID MATCHES FujitsuClang )
    option ( WALBERLA_CXX_COMPILER_IS_FUJITSUCLANG "Use FujitsuClang compiler" ON  )
    include(FujitsuClang)
else()
    option ( WALBERLA_CXX_COMPILER_IS_FUJITSUCLANG "Use FujitsuClang compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_FUJITSUCLANG )

# Check for MPI wrapper
get_filename_component( CXX_COMPILER_WITHOUT_PATH ${CMAKE_CXX_COMPILER} NAME )
if( CXX_COMPILER_WITHOUT_PATH MATCHES "mpi" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "mpi" )
    option ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" ON  )
else()
    option ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER "Compiler is MPI wrapper" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_MPI_WRAPPER )

# Check for intel llvm compiler
if( CMAKE_CXX_COMPILER MATCHES "icpx" OR CMAKE_CXX_COMPILER_ARG1 MATCHES "icpx" )
   option ( WALBERLA_CXX_COMPILER_IS_INTELLLVM "Use Intel LLVM compiler" ON  )
   include(IntelLLVM)
else()
   option ( WALBERLA_CXX_COMPILER_IS_INTELLLVM "Use Intel LLVM compiler" OFF  )
endif()
mark_as_advanced ( WALBERLA_CXX_COMPILER_IS_INTELLLVM )

############################################################################################################################



############################################################################################################################
##
## Compiler Setup
##
############################################################################################################################



# Xcode generator disables -isystem flag, even though current versions of Xcode support it
if(CMAKE_GENERATOR STREQUAL "Xcode")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif()


# disable Xcode 7.3+ linker deduplication pass to speed up linking in debug mode
#if ( APPLE )
#   execute_process( COMMAND ${CMAKE_LINKER} -v OUTPUT_VARIABLE LINKER_VERSION ERROR_VARIABLE LINKER_VERSION )
#   string( REGEX MATCH "ld64-[0-9\\.\\-]+" LINKER_VERSION ${LINKER_VERSION} )
#   string( REGEX MATCHALL "[^\\-]+" LINKER_VERSION ${LINKER_VERSION} )
#   list( GET LINKER_VERSION 0 LINKER_TYPE )
#   list( GET LINKER_VERSION 1 LINKER_VERSION )
#   if( LINKER_TYPE STREQUAL "ld64" AND LINKER_VERSION VERSION_GREATER 264.3.101 )
#       add_flag( CMAKE_EXE_LINKER_FLAGS_DEBUG    "-Wl,-no_deduplicate")
#       add_flag( CMAKE_MODULE_LINKER_FLAGS_DEBUG "-Wl,-no_deduplicate")
#       add_flag( CMAKE_SHARED_LINKER_FLAGS_DEBUG "-Wl,-no_deduplicate")
#   endif()
#endif()


############################################################################################################################


############################################################################################################################
##
## Python
##
#############################################################################################################################

if ( WALBERLA_BUILD_WITH_CODEGEN OR WALBERLA_BUILD_WITH_PYTHON )
   cmake_policy( SET CMP0094 NEW )
   set( Python_FIND_FRAMEWORK LAST )
   find_package( Python COMPONENTS Interpreter Development )
endif ()

############################################################################################################################
##
## Code Generation (pystencils)
##
#############################################################################################################################
if ( WALBERLA_BUILD_WITH_CODEGEN )
   set(LBMPY_MIN_VERSION 1.3.7)
   execute_process(COMMAND ${Python_EXECUTABLE} -c "import lbmpy; print(lbmpy.__version__)"
         RESULT_VARIABLE LBMPY_FOUND OUTPUT_VARIABLE LBMPY_VERSION)
    if(NOT LBMPY_FOUND EQUAL 0)
       message(FATAL_ERROR "WALBERLA_BUILD_WITH_CODEGEN activated but pystencils or lbmpy package not found.
                            Please install lbmpy e.g.: 'pip3 install lbmpy'")
    elseif(LBMPY_VERSION VERSION_LESS LBMPY_MIN_VERSION)
       string(STRIP ${LBMPY_VERSION} LBMPY_VERSION_STRIP)
       message(WARNING
             "lbmpy version ${LBMPY_VERSION_STRIP} was found.\n"
             "We recommend to use at least version ${LBMPY_MIN_VERSION}.")
    endif()
    execute_process(COMMAND ${Python_EXECUTABLE} -c "from pystencils.include import get_pystencils_include_path; print(get_pystencils_include_path())"
                    OUTPUT_VARIABLE PYSTENCILS_INCLUDE_PATH)
    include_directories( ${PYSTENCILS_INCLUDE_PATH} )

    execute_process(COMMAND ${Python_EXECUTABLE} -c "import jinja2"
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

    include(FetchContent)

    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG        v2.13.6
    )

    FetchContent_MakeAvailable(pybind11)

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
## PThread is required in Linux environments by std::thread
##
############################################################################################################################

set( THREADS_PREFER_PTHREAD_FLAG TRUE )
find_package(Threads)
if ( Threads_FOUND )
   if( CMAKE_USE_PTHREADS_INIT )
      add_flag( CMAKE_CXX_FLAGS "-pthread" )
   else()
      add_flag( CMAKE_CXX_FLAGS "${CMAKE_THREAD_LIBS_INIT}" )
   endif()
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
## METIS
##
############################################################################################################################

if ( WALBERLA_BUILD_WITH_PARMETIS )
   # metis is required for parmetis
   set( WALBERLA_BUILD_WITH_METIS TRUE FORCE )
endif ()

if ( WALBERLA_BUILD_WITH_METIS )
   find_package( Metis REQUIRED )

   include_directories( ${METIS_INCLUDE_DIRS} )
   list( APPEND SERVICE_LIBS ${METIS_LIBRARIES} )

endif()

if ( WALBERLA_BUILD_WITH_PARMETIS )
   find_package( Parmetis REQUIRED )

   include_directories( ${PARMETIS_INCLUDE_DIR} )
   list( APPEND SERVICE_LIBS ${PARMETIS_LIBRARY} )
endif()

############################################################################################################################



############################################################################################################################
##
## FFTW3
##
############################################################################################################################

if ( WALBERLA_BUILD_WITH_FFTW )
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
      message("   If OpenMesh required, set OPENMESH_LIBRARY_DIR to the OpenMesh directory.")
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

if ( WALBERLA_BUILD_WITH_OPENMP AND NOT OpenMP_FOUND )
    if( APPLE AND EXISTS /opt/local/lib/libomp AND EXISTS /opt/local/include/libomp ) # find libomp from MacPorts
        set( CMAKE_FRAMEWORK_PATH /opt/local/lib/libomp )
        set( CMAKE_INCLUDE_PATH /opt/local/include/libomp )
    endif()
    find_package( OpenMP )
    if (OpenMP_FOUND)
        add_flag ( CMAKE_C_FLAGS   "${OpenMP_C_FLAGS}" )
        add_flag ( CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}" )
        list ( APPEND SERVICE_LIBS ${OpenMP_CXX_LIBRARIES} )
        if( OpenMP_CXX_INCLUDE_DIRS )
            include_directories( ${OpenMP_CXX_INCLUDE_DIRS} )
        endif()
    else()
         message(FATAL_ERROR "Could NOT enable OpenMP")
    endif()
endif()
############################################################################################################################



############################################################################################################################
##
## CUDA
##
############################################################################################################################
if ( WALBERLA_BUILD_WITH_CUDA )
    if (WALBERLA_BUILD_WITH_HIP)
       message(FATAL_ERROR "For GPU support either use CUDA or HIP. Both simultaneously is not supported.")
    endif()

    include(CheckLanguage)
    check_language(CUDA)
    if( CMAKE_CUDA_COMPILER )

        if(${CMAKE_VERSION} VERSION_GREATER "3.18.0" AND NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES OFF)
            message(WARNING "CMAKE_CUDA_ARCHITECTURES was not set. It is automatically set to: ${CMAKE_CUDA_ARCHITECTURES}")
        endif()

        enable_language(CUDA)

        #include directories and cudart lib is needed for cpp files that use cuda headers/libs
        include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
        find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
        list ( APPEND SERVICE_LIBS ${CUDART_LIBRARY} )

        find_library( NVTX_LIBRARY nvToolsExt ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES} )
        list ( APPEND SERVICE_LIBS ${NVTX_LIBRARY} )

        #CUDA_FOUND is need for our cmake mechanism
        set ( CUDA_FOUND TRUE )
        set (WALBERLA_BUILD_WITH_GPU_SUPPORT TRUE)
    else()
        message( WARNING "CUDA could not be enabled. The host compiler might not be compatible. Check CMakeFiles/CMakeError.log for more information" )
        set ( WALBERLA_BUILD_WITH_CUDA FALSE )
    endif ( )
endif ( )

if (WALBERLA_BUILD_WITH_CUDA AND (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "11.0" OR CMAKE_VERSION VERSION_LESS 3.18.0))
    # CUDA < 11 does not support C++17. std::experimental::any works with C++14, unlike std::any.
    set(CMAKE_CUDA_STANDARD 14)
    set(WALBERLA_USE_STD_EXPERIMENTAL_ANY 1)
endif()

# Can be used in CMake files containing generated files where the file suffix is dependent on the target
if (WALBERLA_BUILD_WITH_CUDA AND CUDA_FOUND)
    set(CODEGEN_FILE_SUFFIX "cu")
else()
    set(CODEGEN_FILE_SUFFIX "cpp")
endif()
############################################################################################################################



############################################################################################################################
##
## ROCm HIP
##
############################################################################################################################
if ( WALBERLA_BUILD_WITH_HIP )
    if (WALBERLA_BUILD_WITH_CUDA)
       message(FATAL_ERROR "For GPU support either use CUDA or HIP. Both simultaneously is not supported.")
    endif()
    if (${CMAKE_VERSION} VERSION_LESS "3.21.0")
       message(FATAL_ERROR "For HIP support CMake > 3.21.0 is needed. Please install a newer version")
    endif()

    include(CheckLanguage)
    check_language(HIP)

    if( CMAKE_HIP_COMPILER )
        enable_language(HIP)
        find_package(hip REQUIRED)
        # since waLBerla also supports CUDA we only use HIP on an AMD platform
        add_compile_definitions(__HIP_PLATFORM_AMD__)
        # include_directories(${HSA_HEADER})
        #HIP_FOUND is need for our cmake mechanism
        set ( HIP_FOUND TRUE )
        set (WALBERLA_BUILD_WITH_GPU_SUPPORT TRUE)
    else()
         message("HIP compiler not found. HIP support is not possible")
        set ( WALBERLA_BUILD_WITH_HIP FALSE )
    endif ( )
endif ( )

if (WALBERLA_BUILD_WITH_HIP AND HIP_FOUND)
    set(CODEGEN_FILE_SUFFIX "hip")
endif()
############################################################################################################################



############################################################################################################################
##
##  Likwid Marker API
##
############################################################################################################################


if ( WALBERLA_BUILD_WITH_LIKWID_MARKERS AND NOT LIKWID_FOUND )
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
    cmake_policy( SET CMP0069 NEW )
    include( CheckIPOSupported )
    check_ipo_supported( RESULT LTO_SUPPORTED LANGUAGES CXX )
    if( LTO_SUPPORTED )
       set( CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE )
    else()
       message( WARNING "Link-time optimization is not supported with this compiler" )
    endif()
endif ( )

############################################################################################################################
##
##  Some more compiler flags that need to happen after any try_compile (e.g. inside FindMPI) (not sure if still true)
##
############################################################################################################################

# Treat warnings as errors
if ( WARNING_ERROR )
   set (CMAKE_COMPILE_WARNING_AS_ERROR ON)
endif ( )

############################################################################################################################
##
##  Half precision
##
############################################################################################################################
if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
   ### Compiler requirements:
   ### Within this project, there are several checks to ensure that the template parameter 'ValueType'
   ### is a floating point number. The check is_floating_point<ValueType> is done primarily in our MPI implementation.
   ### The IEE 754 floating type format _Float16, evaluates to true only if your compiler supports the
   ### open C++23 standard P1467R9 (Extended floating-point types and standard names).
   ### Compare:
   ###  https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p1467r9.html
   ###
   ### Right now (18.12.2023) this is the case only for gcc13.
   ### For more information see:
   ###   https://gcc.gnu.org/projects/cxx-status.html#:~:text=Extended%20floating%2Dpoint%20types%20and%20standard%20names
   ###   https://clang.llvm.org/cxx_status.html#:~:text=Extended%20floating%2Dpoint%20types%20and%20standard%20names

   try_compile( WALBERLA_SUPPORT_HALF_PRECISION "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestFloat16.cpp"
         CXX_STANDARD 23 OUTPUT_VARIABLE TRY_COMPILE_OUTPUT )
   ## message( STATUS ${TRY_COMPILE_OUTPUT} )
   if ( NOT WALBERLA_SUPPORT_HALF_PRECISION )
      message( FATAL_ERROR "Compiler: ${CMAKE_CXX_COMPILER} Version: ${CMAKE_CXX_COMPILER_VERSION} does not support half precision" )
   endif ()

endif () # Check if WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT is set

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
# Subdirectories
#
############################################################################################################################

# Add binary dir for generated headers
include_directories ( ${CMAKE_CURRENT_BINARY_DIR}/src )
# All include paths are specified relative to src/ directory
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/src )

# external
add_subdirectory( extern )

# sources
add_subdirectory ( src )

# Generate file with compile options, and add install rule for it
configure_file ( src/waLBerlaDefinitions.in.h
                 src/waLBerlaDefinitions.h    )

install( FILES ${walberla_BINARY_DIR}/src/waLBerlaDefinitions.h DESTINATION walberla/ )


# test
if ( WALBERLA_BUILD_TESTS )
    add_subdirectory ( tests )
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