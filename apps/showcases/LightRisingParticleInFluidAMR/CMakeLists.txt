add_executable ( LightRisingParticleInFluidAMR LightRisingParticleInFluidAMR.cpp )
target_link_libraries ( LightRisingParticleInFluidAMR walberla::core walberla::mesa_pd walberla::lbm
        walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk
        walberla::geometry walberla::postprocessing )
