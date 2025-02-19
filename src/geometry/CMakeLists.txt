###################################################################################################
#
# Module geometry
#
###################################################################################################

add_library( walberla_geometry )
add_library( walberla::geometry ALIAS walberla_geometry )
target_sources( walberla_geometry
      PRIVATE
      GeometricalFunctions.h
      InitBoundaryHandling.h
      GeometricalFunctions.cpp
      )
target_link_libraries( walberla_geometry PUBLIC walberla::core walberla::domain_decomposition walberla::boundary walberla::field walberla::stencil PRIVATE lodepng )

add_subdirectory( bodies )
add_subdirectory( containment_octree )
add_subdirectory( initializer )
add_subdirectory( mesh )
add_subdirectory( structured )
###################################################################################################
