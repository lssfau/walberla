waLBerla_link_files_to_builddir( *.cfg )

add_executable ( 01_MESA_PD 01_MESAPD.cpp )
target_link_libraries( 01_MESA_PD walberla::blockforest walberla::core walberla::mesa_pd )
