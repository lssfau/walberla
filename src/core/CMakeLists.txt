
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

file ( GLOB_RECURSE sourceFiles ${WALBERLA_GLOB_FILES} )
list( APPEND sourceFiles "${walberla_BINARY_DIR}/src/core/waLBerlaBuildInfo.cpp" )

waLBerla_add_module( FILES ${sourceFiles} EXCLUDE "${walberla_SOURCE_DIR}/src/core/waLBerlaBuildInfo.in.cpp" )
