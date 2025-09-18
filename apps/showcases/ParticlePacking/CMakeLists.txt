waLBerla_link_files_to_builddir( *.cfg )

if (WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
   add_executable ( ParticlePacking ParticlePacking.cpp )
   target_link_libraries( ParticlePacking walberla::blockforest walberla::core walberla::mesa_pd
           walberla::postprocessing walberla::sqlite walberla::vtk )
endif()

