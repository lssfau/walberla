waLBerla_link_files_to_builddir( *.cfg )

if (WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
   waLBerla_add_executable ( NAME ParticlePacking
                             FILES ParticlePacking.cpp
         DEPENDS walberla::blockforest walberla::core walberla::mesa_pd walberla::postprocessing walberla::sqlite walberla::vtk )
endif()

