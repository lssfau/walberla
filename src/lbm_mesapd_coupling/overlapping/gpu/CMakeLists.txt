target_sources(lbm_mesapd_coupling
        PRIVATE
        ParticleAndVolumeFractionMappingKernels.${CODEGEN_FILE_SUFFIX}
        ParticleAndVolumeFractionMappingKernels.h
        )
