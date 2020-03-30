waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable( NAME Mixer
   FILES Mixer.cpp
   DEPENDS core mesa_pd sqlite vtk )
