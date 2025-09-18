waLBerla_link_files_to_builddir( *.cfg )

add_executable( HeatConduction HeatConduction.cpp )
target_link_libraries ( HeatConduction walberla::core walberla::mesa_pd walberla::vtk )
