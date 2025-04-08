###################################################################################################
#
# Tests for LBM_MESAPD_COUPLING functionality
#
###################################################################################################

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_Mapping FILES mapping/ParticleMapping.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_Mapping PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_MovingMapping FILES momentum_exchange_method/MovingParticleMappingMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_MovingMapping PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_DragForceSphere FILES momentum_exchange_method/DragForceSphereMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_DragForceSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_DragForceSphere> --funcTest PROCESSES 2)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs FILES momentum_exchange_method/ForceBetweenTwoStationaryObjectsMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --systemVelocity 0.1 PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --systemVelocity 0.1 PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere FILES momentum_exchange_method/SettlingSphereMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_SettlingSphere> --funcTest PROCESSES 4)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction FILES momentum_exchange_method/PdfReconstructionMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping FILES momentum_exchange_method/UpdateParticleMappingMEM.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field walberla::vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_LubricationCorrection FILES utility/LubricationCorrection.cpp DEPENDS walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_LubricationCorrection PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_VirtualMass FILES utility/VirtualMass.cpp DEPENDS walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_VirtualMass PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_InspectionProbe FILES utility/InspectionProbe.cpp DEPENDS walberla::core walberla::mesa_pd walberla::lbm walberla::lbm_mesapd_coupling walberla::domain_decomposition walberla::field )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_InspectionProbe PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks FILES utility/HydrodynamicForceOnMultipleBlocks.cpp DEPENDS walberla::mesa_pd walberla::lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_EulerAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_VVAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --useVV PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_EulerNoAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --noForceAveraging PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_VVNoAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --noForceAveraging --useVV PROCESSES 4 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping FILES partially_saturated_cells_method/ParticleAndVolumeFractionMappingPSM.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm_mesapd_coupling walberla::stencil walberla::mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping> PROCESSES 27 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping_CPU_GPU FILES partially_saturated_cells_method/codegen/ParticleAndVolumeFractionMapping.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::gpu walberla::field walberla::lbm_mesapd_coupling walberla::stencil walberla::mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping_CPU_GPU COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_ParticleAndVolumeFractionMapping_CPU_GPU> PROCESSES 27 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere FILES partially_saturated_cells_method/DragForceSpherePSM.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::timeloop walberla::mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphereFuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere> --funcTest PROCESSES 8 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere> PROCESSES 8 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere FILES partially_saturated_cells_method/SettlingSpherePSM.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::timeloop walberla::mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphereFuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere> --funcTest PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere> --resolution 70 PROCESSES 4 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphere FILES partially_saturated_cells_method/TorqueSpherePSM.cpp DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::stencil walberla::timeloop )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC1W1FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC1W1   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC1W1SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC1W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC1W1ParallelTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC1W1              PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC2W1FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC2W1   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC2W1SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC2W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC3W1FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC3W1   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC3W1SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC3W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC1W2FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC1W2   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC1W2SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC1W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC2W2FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC2W2   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC2W2SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC2W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC3W2FuncTest     COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --funcTest --SC3W2   PROCESSES 1 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphereSC3W2SingleTest   COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere> --SC3W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )

if (WALBERLA_BUILD_WITH_CODEGEN)
    if (NOT WALBERLA_BUILD_WITH_GPU_SUPPORT OR (WALBERLA_BUILD_WITH_GPU_SUPPORT AND (CMAKE_CUDA_ARCHITECTURES GREATER_EQUAL 60 OR WALBERLA_BUILD_WITH_HIP)))
        waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere_CPU_GPU FILES partially_saturated_cells_method/codegen/DragForceSpherePSM.cpp DEPENDS walberla::blockforest walberla::core walberla::gpu walberla::field walberla::lbm_mesapd_coupling walberla::mesa_pd PSMCodegenPython_srt_sc1 )
        waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere_CPU_GPU_FuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere_CPU_GPU> --funcTest PROCESSES 1)
        waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere_CPU_GPU COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere_CPU_GPU> PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo)

        foreach (collision_setup srt trt-smagorinsky)
            foreach (solid_collision 1 2 3)
                foreach (weighting 1 2)
                    set(config ${collision_setup}_sc${solid_collision}_w${weighting})
                    waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config} FILES partially_saturated_cells_method/codegen/TorqueSpherePSM.cpp DEPENDS walberla::blockforest walberla::core walberla::gpu walberla::field walberla::lbm_mesapd_coupling walberla::mesa_pd PSMCodegenPython_${collision_setup}_sc${solid_collision} )
                    target_compile_definitions(LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config} PRIVATE SC=${solid_collision})
                    target_compile_definitions(LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config} PRIVATE Weighting=${weighting})
                    waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config}_FuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config}> --funcTest PROCESSES 1)
                    waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config} COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_TorqueSphere_CPU_GPU_${config}> PROCESSES 8 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo)
                endforeach ()
            endforeach ()
        endforeach ()

        waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere_CPU_GPU FILES partially_saturated_cells_method/codegen/SettlingSpherePSM.cpp DEPENDS walberla::blockforest walberla::core walberla::gpu walberla::field walberla::lbm walberla::lbm_mesapd_coupling walberla::mesa_pd walberla::timeloop walberla::vtk PSMCodegenPython_trt-smagorinsky_sc2 )
        waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere_CPU_GPU_FuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere_CPU_GPU> --funcTest PROCESSES 1)
        waLBerla_execute_test(NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere_CPU_GPU COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere_CPU_GPU> PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo)
    endif ()
endif ()
