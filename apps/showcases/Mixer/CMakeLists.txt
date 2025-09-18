waLBerla_link_files_to_builddir( *.cfg )

add_executable ( Mixer Mixer.cpp )
target_link_libraries ( Mixer walberla::core walberla::mesa_pd walberla::sqlite walberla::vtk )
