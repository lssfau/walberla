waLBerla_link_files_to_builddir( "*.ipynb" )

add_executable ( IntegratorAccuracy IntegratorAccuracy.cpp )
target_link_libraries ( IntegratorAccuracy walberla::core walberla::mesa_pd )
