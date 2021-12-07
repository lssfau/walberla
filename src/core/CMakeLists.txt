
############################################################################################################################
#
# Build Info
#
############################################################################################################################


include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC WALBERLA_GIT_SHA1)

string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )

set      ( WALBERLA_BUILD_TYPE      ${CMAKE_BUILD_TYPE} )
set      ( WALBERLA_COMPILER_FLAGS  "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${CMAKE_CXX_FLAGS}" )
string   ( REPLACE "\"" "\\\"" WALBERLA_COMPILER_FLAGS ${WALBERLA_COMPILER_FLAGS} )
site_name( WALBERLA_BUILD_MACHINE )

# Generate file with compile options, and add install rule for it
configure_file ( waLBerlaBuildInfo.in.cpp 
                 waLBerlaBuildInfo.cpp    )
                 
############################################################################################################################

include( math/math.cmake )
include( logging/logging.cmake )

add_library( core )
if( MPI_FOUND )
   target_link_libraries( core PUBLIC MPI::MPI_CXX )
endif()
target_link_libraries( core PUBLIC ${SERVICE_LIBS} )
target_sources( core
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
         FunctionTraits.h
         GetPID.h
         Hostname.h
         Macros.h
         MemoryUsage.h
         MultiArrayIO.h
         MultiArrayIO.impl.h
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
add_subdirectory( logging)
add_subdirectory( load_balancing )
add_subdirectory( config )
add_subdirectory( math )
add_subdirectory( timing )
add_subdirectory( grid_generator )