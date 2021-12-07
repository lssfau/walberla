target_sources( mesh
    PRIVATE
    DefaultTesselation.h
    Types.h     
    )
add_subdirectory( communication )
add_subdirectory( raytracing )
add_subdirectory( rigid_body )
add_subdirectory( tesselation )
add_subdirectory( vtk )