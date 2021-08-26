waLBerla_add_executable ( NAME LightRisingParticleInFluidAMR
        FILES LightRisingParticleInFluidAMR.cpp
        DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk geometry postprocessing )
