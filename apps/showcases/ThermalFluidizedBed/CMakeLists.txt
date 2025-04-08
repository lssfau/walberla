waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME ThermalFluidizedBedMEM FILES ThermalFluidizedBedMEM.cpp
                          DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME ThermalFluidizedBedPSM FILES ThermalFluidizedBedPSM.cpp
        DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::gpu walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk MTWCodegenThermalPSM_sc_1 )
target_compile_definitions(ThermalFluidizedBedPSM PRIVATE Weighting=2)
