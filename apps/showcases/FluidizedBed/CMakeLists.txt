waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME FluidizedBedMEM FILES FluidizedBedMEM.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk )

waLBerla_add_executable ( NAME FluidizedBedPSM FILES FluidizedBedPSM.cpp
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk )
