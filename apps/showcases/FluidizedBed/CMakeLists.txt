waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME FluidizedBed FILES FluidizedBed.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd postprocessing timeloop vtk )
