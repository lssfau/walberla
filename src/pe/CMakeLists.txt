
###################################################################################################
#
# Module physics engine (pe)
#
###################################################################################################

add_library( pe )
target_link_libraries( pe PUBLIC core blockforest domain_decomposition geometry stencil vtk PRIVATE lodepng )
target_sources( pe
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