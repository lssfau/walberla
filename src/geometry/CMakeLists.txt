###################################################################################################
#
# Module geometry
#
###################################################################################################

add_library( geometry )
target_sources( geometry
      PRIVATE
      GeometricalFunctions.h
      InitBoundaryHandling.h
      GeometricalFunctions.cpp
      )
target_link_libraries(geometry PUBLIC core domain_decomposition boundary field stencil PRIVATE lodepng)

add_subdirectory( bodies )
add_subdirectory( containment_octree )
add_subdirectory( initializer )
add_subdirectory( mesh )
add_subdirectory( structured )
###################################################################################################
