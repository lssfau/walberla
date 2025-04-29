###################################################################################################
#
# Module mesh
#
###################################################################################################

add_library( walberla_mesh )
add_library( walberla::mesh ALIAS walberla_mesh )
target_link_libraries( walberla_mesh PUBLIC walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::geometry walberla::stencil walberla::mesh_common )
if ( OPENMESH_CORE_FOUND )
   target_link_directories( walberla_mesh PUBLIC OpenMeshCore )
endif()
target_sources( walberla_mesh
      PRIVATE
      MeshConversion.h
      )

add_subdirectory( blockforest )
add_subdirectory( pe )
add_subdirectory( boundary )

###################################################################################################
