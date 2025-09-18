add_library( walberla_mesh_common INTERFACE )
add_library( walberla::mesh_common ALIAS walberla_mesh_common )
target_link_libraries( walberla_mesh_common INTERFACE walberla::core walberla::vtk )

if ( ${OPENMESH_FOUND} )
    target_link_libraries ( walberla_mesh_common INTERFACE ${OPENMESH_LIBRARIES} )
endif()

target_sources( walberla_mesh_common
      INTERFACE
      DistanceFunction.h
      QHull.h
      MeshBodyWallDistanceCallback.h
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