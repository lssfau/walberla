###################################################################################################
#
# Module mesh
#
###################################################################################################

add_library( mesh )
target_link_libraries( mesh PUBLIC blockforest boundary core domain_decomposition field geometry pe stencil mesh_common )
if ( OPENMESH_CORE_FOUND )
   target_link_directories( mesh PUBLIC OpenMeshCore)
endif()
target_sources( mesh
      PRIVATE
      MeshConversion.h
      )

add_subdirectory( blockforest )
add_subdirectory( pe )
add_subdirectory( boundary )

###################################################################################################
