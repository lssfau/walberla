waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable( NAME HeatConduction
   FILES HeatConduction.cpp
      DEPENDS walberla::core walberla::mesa_pd walberla::vtk )
