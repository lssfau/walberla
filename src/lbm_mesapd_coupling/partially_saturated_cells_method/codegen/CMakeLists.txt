target_sources(lbm_mesapd_coupling
        PRIVATE
        ParticleAndVolumeFractionMappingSweepsCPU.h
        PSMSweepCollection.h
        PSMWrapperSweepsCPU.h
)
if (WALBERLA_BUILD_WITH_CODEGEN)
    foreach (collision_setup srt trt mrt cumulant srt-smagorinsky trt-smagorinsky)
        foreach (solid_collision 1 2 3)
            set(config ${collision_setup}_sc${solid_collision})
            walberla_generate_target_from_python(NAME PSMCodegenPython_${config}
                    FILE PSMCodegen.py
                    CODEGEN_CFG ${config}_${MAX_PARTICLES_PER_CELL}
                    OUT_FILES LBMSweep.${CODEGEN_FILE_SUFFIX} LBMSweep.h LBMSplitSweep.${CODEGEN_FILE_SUFFIX} LBMSplitSweep.h
                    PSMSweep.${CODEGEN_FILE_SUFFIX} PSMSweep.h PSMSweepSplit.${CODEGEN_FILE_SUFFIX} PSMSweepSplit.h
                    PSMPackInfo.${CODEGEN_FILE_SUFFIX} PSMPackInfo.h InitializeDomainForPSM.${CODEGEN_FILE_SUFFIX} InitializeDomainForPSM.h
                    PSM_NoSlip.${CODEGEN_FILE_SUFFIX} PSM_NoSlip.h PSM_UBB.${CODEGEN_FILE_SUFFIX} PSM_UBB.h PSM_Density.${CODEGEN_FILE_SUFFIX} PSM_Density.h
                    PSM_FreeSlip.${CODEGEN_FILE_SUFFIX} PSM_FreeSlip.h PSM_InfoHeader.h PSM_MacroGetter.cpp PSM_MacroGetter.h
                    PSM_MacroSetter.cpp PSM_MacroSetter.h)
            add_dependencies(PSMCodegenPython_${config} MAX_PARTICLES_PER_CELL)
        endforeach ()
    endforeach ()
endif ()
if (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP))
    target_sources(lbm_mesapd_coupling
            PRIVATE
            ParticleAndVolumeFractionMappingSweepsGPU.h
            ParticleAndVolumeFractionMappingKernels.${CODEGEN_FILE_SUFFIX}
            ParticleAndVolumeFractionMappingKernels.h
            PSMUtilityGPU.h
            PSMWrapperKernels.${CODEGEN_FILE_SUFFIX}
            PSMWrapperKernels.h
            PSMWrapperSweepsGPU.h
    )
endif ()
