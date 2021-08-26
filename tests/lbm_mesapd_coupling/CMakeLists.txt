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

