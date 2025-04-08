waLBerla_add_executable ( NAME LightRisingParticleInFluidAMR
        FILES LightRisingParticleInFluidAMR.cpp
      DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk walberla::geometry walberla::postprocessing )
