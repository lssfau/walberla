waLBerla_link_files_to_builddir( *.cfg )

if (OPENMESH_FOUND AND WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
    waLBerla_add_executable ( NAME MESA_PD_MESH_APP_PegIntoSphereBed
            FILES PegIntoSphereBed.cpp
            DEPENDS blockforest mesh_common mesa_pd core vtk )
endif()