waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable ( NAME 01_MESA_PD
   FILES 01_MESAPD.cpp
                          DEPENDS blockforest core mesa_pd )
