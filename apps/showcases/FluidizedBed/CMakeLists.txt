waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME FluidizedBedMEM FILES FluidizedBedMEM.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd timeloop vtk )

waLBerla_add_executable ( NAME FluidizedBedPSM FILES FluidizedBedPSM.cpp
        DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling mesa_pd timeloop vtk )
