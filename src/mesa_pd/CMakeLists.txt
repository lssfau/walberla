
###################################################################################################
#
# Module rigid particle dynamics (RPD)
#
###################################################################################################

add_library( mesa_pd )
target_link_libraries( mesa_pd PUBLIC blockforest core stencil vtk )
add_subdirectory( collision_detection )
add_subdirectory( common )
add_subdirectory( data )
add_subdirectory( domain )
add_subdirectory( kernel )
add_subdirectory( mpi )
add_subdirectory( sorting )
add_subdirectory( vtk )
if(OPENMESH_CORE_FOUND)
    set( WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE ON CACHE INTERNAL "")
    target_link_libraries( mesa_pd PUBLIC mesh_common )
    message( STATUS "MESA-PD: ConvexPolyhedron shape is available (OpenMesh dependency satisfied)" )
else()
    set( WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE OFF CACHE INTERNAL "")
    message( STATUS "MESA-PD: ConvexPolyhedron shape is unavailable (OpenMesh not found)" )
endif()