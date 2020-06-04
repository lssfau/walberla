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

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects FILES momentum_exchange_method/ForceBetweenTwoStationaryObjects.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSS1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSS2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSS3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSS4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --systemVelocity 0.1 PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSW1 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useSphereWallSetup PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSW2 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useSphereWallSetup --useSBB PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSW3 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useSphereWallSetup --useCompressible PROCESSES 1)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjectsSW4 COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_ForceBetweenTwoStationaryObjects> --useSphereWallSetup --systemVelocity 0.1 PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere FILES momentum_exchange_method/SettlingSphere.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_SettlingSphere COMMAND $<TARGET_FILE:LBM_MESAPD_COUPLING_MEM_SettlingSphere> --funcTest PROCESSES 4)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction FILES momentum_exchange_method/PdfReconstruction.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_PdfReconstruction PROCESSES 3)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping FILES momentum_exchange_method/UpdateParticleMapping.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field vtk)
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_MEM_UpdateParticleMapping PROCESSES 1)

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTILITY_LubricationCorrection FILES utility/LubricationCorrection.cpp DEPENDS mesa_pd lbm_mesapd_coupling )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTILITY_LubricationCorrection PROCESSES 1 )

waLBerla_compile_test( NAME LBM_MESAPD_COUPLING_UTILITY_InspectionProbe FILES utility/InspectionProbe.cpp DEPENDS core mesa_pd lbm lbm_mesapd_coupling domain_decomposition field )
waLBerla_execute_test( NAME LBM_MESAPD_COUPLING_UTILITY_InspectionProbe PROCESSES 1 )

