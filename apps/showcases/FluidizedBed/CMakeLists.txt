waLBerla_link_files_to_builddir( "*.prm" )

add_executable ( FluidizedBedMEM FluidizedBedMEM.cpp )
target_link_libraries ( FluidizedBedMEM walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd
        walberla::timeloop walberla::vtk )

add_executable ( FluidizedBedPSM FluidizedBedPSM.cpp )
target_link_libraries ( FluidizedBedPSM walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition
        walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk )
