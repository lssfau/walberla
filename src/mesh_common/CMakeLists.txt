###################################################################################################
#
# Module mesh_common
#
###################################################################################################

add_library( mesh_common INTERFACE )
target_link_libraries( mesh_common INTERFACE core vtk )
if ( OPENMESH_CORE_FOUND )
   target_link_directories( mesh_common INTERFACE OpenMeshCore)
endif()
target_sources( mesh_common
      INTERFACE
      DistanceFunction.h
      QHull.h
      MeshOperations.h
      MeshIO.h
      DistanceComputations.h
      TriangleMeshes.h
      MatrixVectorOperations.h
      OpenMeshBufferTypeExtensions.h
      PolyMeshes.h
      )

add_subdirectory( vtk )
add_subdirectory( distance_octree )

###################################################################################################
