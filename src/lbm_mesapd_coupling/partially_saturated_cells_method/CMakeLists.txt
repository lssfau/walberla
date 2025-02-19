target_sources( walberla_lbm_mesapd_coupling
    PRIVATE
    PSMSweep.h
    ParticleAndVolumeFractionMapping.h
    PSMUtility.h
    )

add_subdirectory(codegen)
