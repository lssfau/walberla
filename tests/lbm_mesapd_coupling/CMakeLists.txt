###################################################################################################
#
# Tests for LBM_MESAPD_COUPLING functionality
#
###################################################################################################

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_Mapping FILES mapping/ParticleMapping.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_Mapping PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_MovingMapping FILES momentum_exchange_method/MovingParticleMapping.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_MovingMapping PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_DragForceSphere FILES momentum_exchange_method/DragForceSphere.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_DragForceSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_DragForceSphere> --funcTest PROCESSES 2)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs FILES momentum_exchange_method/ForceBetweenTwoStationaryObjects.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSS4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --systemVelocity 0.1 PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjsSW4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceTwoStatObjs> --useSphereWallSetup --systemVelocity 0.1 PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere FILES momentum_exchange_method/SettlingSphere.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_SettlingSphere> --funcTest PROCESSES 4)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction FILES momentum_exchange_method/PdfReconstruction.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping FILES momentum_exchange_method/UpdateParticleMapping.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_LubricationCorrection FILES utility/LubricationCorrection.cpp DEPENDS mesa_pd lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_LubricationCorrection PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_VirtualMass FILES utility/VirtualMass.cpp DEPENDS mesa_pd lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_VirtualMass PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_InspectionProbe FILES utility/InspectionProbe.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_InspectionProbe PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks FILES utility/HydrodynamicForceOnMultipleBlocks.cpp DEPENDS mesa_pd lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_EulerAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_VVAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --useVV PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_EulerNoAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --noForceAveraging PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks_VVNoAvg COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_UTIL_HydForceMultBlocks> --noForceAveraging --useVV PROCESSES 4 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_BodyAndVolumeFractionMapping FILES partially_saturated_cells_method/ParticleAndVolumeFractionMapping.cpp DEPENDS blockforest boundary core field lbm_mesapd_coupling stencil mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_BodyAndVolumeFractionMapping PROCESSES 27 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere FILES partially_saturated_cells_method/DragForceSphere.cpp DEPENDS blockforest boundary core field lbm lbm_mesapd_coupling timeloop mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphereFuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere> --funcTest PROCESSES 8 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_DragForceSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_DragForceSphere> PROCESSES 8 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere FILES partially_saturated_cells_method/SettlingSphere.cpp DEPENDS blockforest boundary core domain_decomposition field lbm lbm_mesapd_coupling timeloop mesa_pd )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphereFuncTest COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere> --funcTest PROCESSES 4 )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_PSM_SettlingSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_PSM_SettlingSphere> --resolution 70 PROCESSES 4 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_PSM_TorqueSphere FILES partially_saturated_cells_method/TorqueSphere.cpp DEPENDS blockforest boundary core domain_decomposition field lbm stencil timeloop )
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
