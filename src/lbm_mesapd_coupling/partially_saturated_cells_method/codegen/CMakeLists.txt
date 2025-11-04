target_sources( walberla_lbm_mesapd_coupling
        PRIVATE
        ParticleAndVolumeFractionMappingSweepsCPU.h
        PSMSweepCollection.h
        PSMWrapperSweepsCPU.h
)
if (WALBERLA_BUILD_WITH_CODEGEN)
    foreach (collision_setup srt trt mrt cumulant srt-smagorinsky trt-smagorinsky)
        foreach (solid_collision 1 2 3)
            foreach (force 0 1)
                set(config ${collision_setup}_sc${solid_collision})
                if (force EQUAL 1)
                    set(config "${config}_force")
                endif ()
                walberla_generate_target_from_python(NAME PSMCodegenPython_${config}
                        FILE PSMCodegen.py
                        CODEGEN_CFG ${config}_${MAX_PARTICLES_PER_CELL}
                        OUT_FILES PSMSweep.${CODEGEN_FILE_SUFFIX} PSMSweep.h
                        PSM_MacroSetter.${CODEGEN_FILE_SUFFIX} PSM_MacroSetter.h PSM_MacroGetter.h PSM_MacroGetter.${CODEGEN_FILE_SUFFIX}
                        LBM_MacroSetter.${CODEGEN_FILE_SUFFIX} LBM_MacroSetter.h
                        LBMSweepCollection.h LBMSweepCollection.${CODEGEN_FILE_SUFFIX}
                        LBMStorageSpecification.h LBMStorageSpecification.${CODEGEN_FILE_SUFFIX}
                        LBMBoundaryCollection.h
                        LBM_UBB.h LBM_UBB.${CODEGEN_FILE_SUFFIX}
                        LBM_FreeSlip.h LBM_FreeSlip.${CODEGEN_FILE_SUFFIX}
                        LBM_NoSlip.h LBM_NoSlip.${CODEGEN_FILE_SUFFIX}
                        LBM_FixedDensity.h LBM_FixedDensity.${CODEGEN_FILE_SUFFIX}
                        LBM_FixedDensity_2.h LBM_FixedDensity_2.${CODEGEN_FILE_SUFFIX}
                        LBM_InfoHeader.h
                )
                add_dependencies(PSMCodegenPython_${config} MAX_PARTICLES_PER_CELL)
            endforeach ()
        endforeach ()
    endforeach ()
endif ()
if (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP))
    target_sources( walberla_lbm_mesapd_coupling
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
