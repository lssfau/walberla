waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable( NAME HeatConduction
   FILES HeatConduction.cpp
   DEPENDS core mesa_pd vtk )
