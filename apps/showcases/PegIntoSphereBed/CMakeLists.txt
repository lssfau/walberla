waLBerla_link_files_to_builddir( *.cfg )

if (OPENMESH_FOUND AND WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
    add_executable ( MESA_PD_MESH_APP_PegIntoSphereBed PegIntoSphereBed.cpp )
    target_link_libraries ( MESA_PD_MESH_APP_PegIntoSphereBed walberla::mesh_common walberla::mesa_pd walberla::core
            walberla::vtk )
endif()