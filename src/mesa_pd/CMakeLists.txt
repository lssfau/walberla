
###################################################################################################
#
# Module rigid particle dynamics (RPD)
#
###################################################################################################

add_library( walberla_mesa_pd )
add_library( walberla::mesa_pd ALIAS walberla_mesa_pd )
target_link_libraries( walberla_mesa_pd PUBLIC walberla::blockforest walberla::core walberla::stencil walberla::vtk )
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
    target_link_libraries( walberla_mesa_pd PUBLIC walberla::mesh_common )
    message( DEBUG "MESA-PD: ConvexPolyhedron shape is available (OpenMesh dependency satisfied)" )
else()
    set( WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE OFF CACHE INTERNAL "")
    message( DEBUG "MESA-PD: ConvexPolyhedron shape is unavailable (OpenMesh not found)" )
endif()