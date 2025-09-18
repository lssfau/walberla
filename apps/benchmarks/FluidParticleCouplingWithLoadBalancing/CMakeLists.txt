add_executable ( FluidParticleWorkloadEvaluation FluidParticleWorkloadEvaluation.cpp )
target_link_libraries ( FluidParticleWorkloadEvaluation walberla::blockforest walberla::boundary walberla::core
        walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing
        walberla::timeloop walberla::vtk )

add_executable( FluidParticleWorkloadDistribution FluidParticleWorkloadDistribution.cpp )
target_link_libraries ( FluidParticleWorkloadDistribution walberla::blockforest walberla::boundary walberla::core
        walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::postprocessing
        walberla::timeloop walberla::vtk )
