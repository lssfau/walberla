
############################################################################################################################
#
# Build Info
#
############################################################################################################################


include( GetGitRevisionDescription )
get_git_head_revision( GIT_REFSPEC WALBERLA_GIT_SHA1 )

string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )

set( WALBERLA_BUILD_TYPE ${CMAKE_BUILD_TYPE} )
set( WALBERLA_COMPILER_FLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${CMAKE_CXX_FLAGS}" )
string( REPLACE "\"" "\\\"" WALBERLA_COMPILER_FLAGS ${WALBERLA_COMPILER_FLAGS} )
site_name( WALBERLA_BUILD_MACHINE )

# Generate file with compile options, and add install rule for it
configure_file( waLBerlaBuildInfo.in.cpp
      waLBerlaBuildInfo.cpp )

############################################################################################################################

include( math/math.cmake )
include( logging/logging.cmake )

add_library( walberla_core )
add_library( walberla::core ALIAS walberla_core )
if ( MPI_FOUND )
   target_link_libraries( walberla_core PUBLIC MPI::MPI_CXX )
endif ()

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
   # Actual support for float16 is available only since C++23
   #   before is_arithmetic and is_floating_point evaluated to false,
   #   also many STL functions are compatible with float16 only since C++23.
   # Which features are actually supported depend on the compiler
   target_compile_features( walberla_core PUBLIC cxx_std_23 )
endif ()

target_link_libraries( walberla_core PUBLIC ${SERVICE_LIBS} )
target_include_directories( walberla_core PUBLIC ${walberla_SOURCE_DIR}/src/ ${walberla_BINARY_DIR}/src/ )
target_sources( walberla_core
      PRIVATE
      Abort.cpp
      Conversion.cpp
      DataTypes.cpp
      Environment.cpp
      GetPID.cpp
      Hostname.cpp
      MemoryUsage.cpp
      RandomUUID.cpp
      Sleep.cpp
      waLBerlaBuildInfo.cpp
      Abort.h
      AllSet.h
      Any.h
      Array.h
      ConcatIterator.h
      Conversion.h
      DataTypes.h
      EndianIndependentSerialization.h
      Environment.h
      Filesystem.h
      Format.hpp
      FunctionTraits.h
      GetPID.h
      Hostname.h
      Macros.h
      MemoryUsage.h
      NonCopyable.h
      NonCreateable.h
      OpenMP.h
      Optional.h
      RandomUUID.h
      Regex.h
      Sanitizer.h
      Set.h
      SharedFunctor.h
      Sleep.h
      STLIO.h
      stringToNum.h
      StringUtility.h
      StringUtility.impl.h
      typeToString.h
      UniqueID.h
      Variant.h
      VectorTrait.h
      waLBerlaBuildInfo.h
)

add_subdirectory( cell )
add_subdirectory( debug )
add_subdirectory( mpi )
add_subdirectory( logging )
add_subdirectory( load_balancing )
add_subdirectory( config )
add_subdirectory( math )
add_subdirectory( timing )
add_subdirectory( grid_generator )