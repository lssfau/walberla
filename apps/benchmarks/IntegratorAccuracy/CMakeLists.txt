waLBerla_link_files_to_builddir( "*.ipynb" )

waLBerla_add_executable ( NAME IntegratorAccuracy 
                          DEPENDS core mesa_pd )
