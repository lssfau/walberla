add_library( walberla_lbm_generated )
add_library( walberla::lbm_generated ALIAS walberla_lbm_generated )

target_link_libraries( walberla_lbm_generated
        PUBLIC
      walberla::blockforest
      walberla::boundary
      walberla::communication
      walberla::core
      walberla::domain_decomposition
      walberla::field
      walberla::geometry
      walberla::stencil
      walberla::timeloop
      walberla::vtk
        )

add_subdirectory( boundary )
add_subdirectory( communication )
add_subdirectory( gpu )
add_subdirectory( evaluation )
add_subdirectory( field )
add_subdirectory( refinement )
add_subdirectory( storage_specification )
add_subdirectory( sweep_collection )