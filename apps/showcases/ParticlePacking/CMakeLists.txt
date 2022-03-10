waLBerla_link_files_to_builddir( *.cfg )

if (WALBERLA_MESAPD_CONVEX_POLYHEDRON_AVAILABLE)
   waLBerla_add_executable ( NAME ParticlePacking
                             FILES ParticlePacking.cpp
                             DEPENDS blockforest core mesa_pd postprocessing sqlite vtk )
endif()

