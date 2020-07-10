
###################################################################################################
#
# Module rigid particle dynamics (RPD)
#
###################################################################################################

waLBerla_add_module( DEPENDS blockforest core stencil vtk OPTIONAL_DEPENDS mesh_common )
if(OPENMESH_CORE_FOUND)
    set( WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE ON CACHE INTERNAL "")
    message( STATUS "MESA-PD: ConvexPolyhedron shape is available (OpenMesh dependency satisfied)" )
else()
    set( WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE OFF CACHE INTERNAL "")
    message( STATUS "MESA-PD: ConvexPolyhedron shape is unavailable (OpenMesh not found)" )
endif()