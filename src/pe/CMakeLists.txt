
###################################################################################################
#
# Module physics engine (pe)
#
###################################################################################################

add_library( walberla_pe )
add_library( walberla::pe ALIAS walberla_pe )
target_link_libraries( walberla_pe PUBLIC walberla::core walberla::blockforest walberla::domain_decomposition walberla::geometry walberla::stencil walberla::vtk PRIVATE lodepng )
target_sources( walberla_pe
      PRIVATE
      Config.h
      Thresholds.h
      Thresholds.cpp
      Materials.cpp
      Types.h
      basic.h
      BlockFunctions.h
      Materials.h
      )

add_subdirectory( cr )
add_subdirectory( amr )
add_subdirectory( vtk )
add_subdirectory( raytracing )
add_subdirectory( synchronization )
add_subdirectory( ccd )
add_subdirectory( debug )
add_subdirectory( communication )
add_subdirectory( bg )
add_subdirectory( utility )
add_subdirectory( contact )
add_subdirectory( statistics )
add_subdirectory( collision )
add_subdirectory( rigidbody )
add_subdirectory( fcd )