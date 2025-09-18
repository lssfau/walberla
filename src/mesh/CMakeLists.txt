add_library( walberla_mesh )
add_library( walberla::mesh ALIAS walberla_mesh )
target_link_libraries( walberla_mesh PUBLIC walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::geometry walberla::stencil walberla::mesh_common )

if ( ${OPENMESH_FOUND} )
    target_link_directories( walberla_mesh PUBLIC ${OPENMESH_LIBRARIES} )
endif ()

target_sources( walberla_mesh
        PUBLIC
        MeshConversion.h
)

add_subdirectory( blockforest )
add_subdirectory( boundary )