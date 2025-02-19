waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable( NAME Mixer
   FILES Mixer.cpp
      DEPENDS walberla::core walberla::mesa_pd walberla::sqlite walberla::vtk )
