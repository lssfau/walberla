waLBerla_link_files_to_builddir( *.cfg )

waLBerla_add_executable ( NAME 01_LennardJones
                          FILES 01_LennardJones.cpp
                          DEPENDS core mesa_pd )
