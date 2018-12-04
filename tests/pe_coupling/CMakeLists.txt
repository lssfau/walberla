###################################################################################################
#
# Tests for pe_coupling module
#
###################################################################################################


###################################################################################################
# Partially saturated cells method tests
###################################################################################################

waLBerla_compile_test( FILES partially_saturated_cells_method/DragForceSpherePSM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME DragForceSpherePSMSC1W1FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC1W1 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC1W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC1W1            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC1W1ParallelTest COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC1W1            PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC2W1FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC2W1 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC2W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC2W1            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC3W1FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC3W1 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC3W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC3W1            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC1W2FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC1W2 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC1W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC1W2            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC2W2FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC2W2 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC2W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC2W2            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMSC3W2FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC3W2 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMSC3W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSM> --PSMVariant SC3W2            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES partially_saturated_cells_method/DragForceSpherePSMRefinement.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC1W1FuncTest     COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC1W1 --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC1W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC1W1            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC1W1ParallelTest COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC1W1            PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC2W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC2W1            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC3W1SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC3W1            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC1W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC1W2            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC2W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC2W2            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSpherePSMRefinementSC3W2SingleTest   COMMAND $<TARGET_FILE:DragForceSpherePSMRefinement> --PSMVariant SC3W2            PROCESSES 1 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES partially_saturated_cells_method/SegreSilberbergPSM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME SegreSilberbergPSMSC1W1FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC1W1 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC1W1Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC1W1             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergPSMSC2W1FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC2W1 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC2W1Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC2W1             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergPSMSC3W1FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC3W1 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC3W1Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC3W1             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergPSMSC1W2FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC1W2 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC1W2Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC1W2             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergPSMSC2W2FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC2W2 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC2W2Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC2W2             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergPSMSC3W2FuncTest     COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC3W2 --funcTest  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergPSMSC3W2Test         COMMAND $<TARGET_FILE:SegreSilberbergPSM> --SC3W2             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES partially_saturated_cells_method/TorqueSpherePSM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME TorqueSpherePSMSC1W1FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC1W1   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC1W1SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC1W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC1W1ParallelTest COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC1W1              PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC2W1FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC2W1   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC2W1SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC2W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC3W1FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC3W1   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC3W1SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC3W1              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC1W2FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC1W2   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC1W2SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC1W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC2W2FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC2W2   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC2W2SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC2W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSpherePSMSC3W2FuncTest     COMMAND $<TARGET_FILE:TorqueSpherePSM> --funcTest --SC3W2   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSpherePSMSC3W2SingleTest   COMMAND $<TARGET_FILE:TorqueSpherePSM> --SC3W2              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )

###################################################################################################
# Momentum exchange method tests
###################################################################################################

waLBerla_compile_test( FILES momentum_exchange_method/BodyAtBlockBoarderCheck.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME BodyAtBlockBoarderCheckTest1        COMMAND $<TARGET_FILE:BodyAtBlockBoarderCheck> 1         PROCESSES 8 )
waLBerla_execute_test( NAME BodyAtBlockBoarderCheckTest2        COMMAND $<TARGET_FILE:BodyAtBlockBoarderCheck> 2         PROCESSES 8 LABELS longrun )
waLBerla_execute_test( NAME BodyAtBlockBoarderCheckTest3        COMMAND $<TARGET_FILE:BodyAtBlockBoarderCheck> 3         PROCESSES 8 LABELS verylongrun )

waLBerla_compile_test( FILES momentum_exchange_method/BodyMappingTest.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME BodyMappingTest COMMAND $<TARGET_FILE:BodyMappingTest> PROCESSES 1 )

waLBerla_compile_test( FILES momentum_exchange_method/DragForceSphereMEM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME DragForceSphereMEMFuncTest        COMMAND $<TARGET_FILE:DragForceSphereMEM> --funcTest          PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSphereMEMSingleTest      COMMAND $<TARGET_FILE:DragForceSphereMEM>                     PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMParallelTest    COMMAND $<TARGET_FILE:DragForceSphereMEM>                     PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMCLIFuncTest     COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant CLI --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSphereMEMCLISingleTest   COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant CLI            PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMCLIParallelTest COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant CLI            PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMMRFuncTest      COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant MR  --funcTest PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSphereMEMMRSingleTest    COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant MR             PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMMRParallelTest  COMMAND $<TARGET_FILE:DragForceSphereMEM> --MEMVariant MR             PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES momentum_exchange_method/DragForceSphereMEMRefinement.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME DragForceSphereMEMRefinementFuncTest        COMMAND $<TARGET_FILE:DragForceSphereMEMRefinement> --funcTest          PROCESSES 1 )
waLBerla_execute_test( NAME DragForceSphereMEMRefinementSingleTest      COMMAND $<TARGET_FILE:DragForceSphereMEMRefinement>                     PROCESSES 4 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMRefinementParallelTest    COMMAND $<TARGET_FILE:DragForceSphereMEMRefinement>                     PROCESSES 5 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME DragForceSphereMEMRefinementCLITest         COMMAND $<TARGET_FILE:DragForceSphereMEMRefinement> --MO_CLI            PROCESSES 4 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES momentum_exchange_method/GlobalBodyAsBoundaryMEMStaticRefinement.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME GlobalBodyAsBoundaryMEMStaticRefinementTest        COMMAND $<TARGET_FILE:GlobalBodyAsBoundaryMEMStaticRefinement>          PROCESSES 1 )

waLBerla_compile_test( FILES momentum_exchange_method/LubricationCorrectionMEM.cpp DEPENDS blockforest timeloop )
waLBerla_execute_test( NAME LubricationCorrectionMEMFuncTest         COMMAND $<TARGET_FILE:LubricationCorrectionMEM> --funcTest                   PROCESSES 3 )
waLBerla_execute_test( NAME LubricationCorrectionMEMSphereSphereTest COMMAND $<TARGET_FILE:LubricationCorrectionMEM> --split --fzyx --sphSphTest  PROCESSES 6 LABELS longrun             CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME LubricationCorrectionMEMSphereWallTest   COMMAND $<TARGET_FILE:LubricationCorrectionMEM> --split --fzyx --sphWallTest PROCESSES 3 LABELS longrun verylongrun CONFIGURATIONS Release RelWithDbgInfo )
  
waLBerla_compile_test( FILES momentum_exchange_method/PeriodicParticleChannelMEM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME PeriodicParticleChannelMEMTest COMMAND $<TARGET_FILE:PeriodicParticleChannelMEM> --shortrun PROCESSES 4 LABELS longrun )

waLBerla_compile_test( FILES momentum_exchange_method/TaylorCouetteFlowMEM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME TaylorCouetteFlowMEMTest COMMAND $<TARGET_FILE:TaylorCouetteFlowMEM> --shortrun )
waLBerla_execute_test( NAME TaylorCouetteFlowMEMTest2 COMMAND $<TARGET_FILE:TaylorCouetteFlowMEM>  PROCESSES 4 LABELS verylongrun)

waLBerla_compile_test( FILES momentum_exchange_method/SegreSilberbergMEM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME SegreSilberbergMEMSyncNextNeighbor    COMMAND $<TARGET_FILE:SegreSilberbergMEM> --shortrun                    PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMSyncShadowOwner     COMMAND $<TARGET_FILE:SegreSilberbergMEM> --shortrun --syncShadowOwners PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMTest                COMMAND $<TARGET_FILE:SegreSilberbergMEM>                      PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergMEMCLIFuncTest         COMMAND $<TARGET_FILE:SegreSilberbergMEM> --funcTest --MO_CLI  PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMCLITest             COMMAND $<TARGET_FILE:SegreSilberbergMEM> --MO_CLI             PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergMEMMRFuncTest          COMMAND $<TARGET_FILE:SegreSilberbergMEM> --funcTest --MO_MR   PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMMRTest              COMMAND $<TARGET_FILE:SegreSilberbergMEM> --MO_MR              PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergMEMEanReconFuncTest    COMMAND $<TARGET_FILE:SegreSilberbergMEM> --funcTest --eanReconstructor PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMEanReconTest        COMMAND $<TARGET_FILE:SegreSilberbergMEM> --MO_CLI --eanReconstructor   PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SegreSilberbergMEMExtReconFuncTest    COMMAND $<TARGET_FILE:SegreSilberbergMEM> --funcTest --extReconstructor PROCESSES 9 )
waLBerla_execute_test( NAME SegreSilberbergMEMExtReconTest        COMMAND $<TARGET_FILE:SegreSilberbergMEM> --MO_CLI --extReconstructor   PROCESSES 18 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES momentum_exchange_method/SettlingSphereMEM.cpp DEPENDS timeloop )
waLBerla_execute_test( NAME SettlingSphereMEMFuncTestSerial      COMMAND $<TARGET_FILE:SettlingSphereMEM> --funcTest PROCESSES 1 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME SettlingSphereMEMFuncTestParallel    COMMAND $<TARGET_FILE:SettlingSphereMEM> --funcTest PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_compile_test( FILES momentum_exchange_method/SettlingSphereMEMDynamicRefinement.cpp DEPENDS timeloop )
waLBerla_execute_test( NAME SettlingSphereMEMDynamicRefinementFuncTestSerial      COMMAND $<TARGET_FILE:SettlingSphereMEMDynamicRefinement> --funcTest PROCESSES 1 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test( NAME SettlingSphereMEMDynamicRefinementFuncTestParallel    COMMAND $<TARGET_FILE:SettlingSphereMEMDynamicRefinement> --funcTest PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_compile_test( FILES momentum_exchange_method/SettlingSphereMEMStaticRefinement.cpp DEPENDS timeloop )
waLBerla_execute_test( NAME SettlingSphereMEMStaticRefinementFuncTestSerial      COMMAND $<TARGET_FILE:SettlingSphereMEMStaticRefinement> --funcTest PROCESSES 1 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test( NAME SettlingSphereMEMStaticRefinementFuncTestParallel    COMMAND $<TARGET_FILE:SettlingSphereMEMStaticRefinement> --funcTest PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_compile_test( FILES momentum_exchange_method/TorqueSphereMEM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME TorqueSphereMEMFuncTest        COMMAND $<TARGET_FILE:TorqueSphereMEM> --funcTest           PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSphereMEMSingleTest      COMMAND $<TARGET_FILE:TorqueSphereMEM>                      PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSphereMEMParallelTest    COMMAND $<TARGET_FILE:TorqueSphereMEM>                      PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSphereMEMCLIFuncTest     COMMAND $<TARGET_FILE:TorqueSphereMEM> --funcTest --MO_CLI  PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSphereMEMCLISingleTest   COMMAND $<TARGET_FILE:TorqueSphereMEM> --MO_CLI             PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSphereMEMCLIParallelTest COMMAND $<TARGET_FILE:TorqueSphereMEM> --MO_CLI             PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSphereMEMMRFuncTest      COMMAND $<TARGET_FILE:TorqueSphereMEM> --funcTest --MO_MR   PROCESSES 1 )
waLBerla_execute_test( NAME TorqueSphereMEMMRSingleTest    COMMAND $<TARGET_FILE:TorqueSphereMEM> --MO_MR              PROCESSES 1 LABELS longrun     CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME TorqueSphereMEMMRParallelTest  COMMAND $<TARGET_FILE:TorqueSphereMEM> --MO_MR              PROCESSES 8 LABELS verylongrun CONFIGURATIONS Release RelWithDbgInfo )  

waLBerla_compile_test( FILES momentum_exchange_method/SquirmerTest.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME SquirmerShortTest COMMAND $<TARGET_FILE:SquirmerTest> --shortrun PROCESSES 1 )
waLBerla_execute_test( NAME SquirmerTest COMMAND $<TARGET_FILE:SquirmerTest> PROCESSES 1 CONFIGURATIONS Release RelWithDbgInfo  )

###################################################################################################
# Discrete particle methods tests
###################################################################################################

waLBerla_compile_test( FILES discrete_particle_methods/SphereWallCollisionBehaviorDPM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME SphereWallCollisionBehaviorDPMFuncTest COMMAND $<TARGET_FILE:SphereWallCollisionBehaviorDPM> --funcTest PROCESSES 1 )

waLBerla_compile_test( FILES discrete_particle_methods/HinderedSettlingDynamicsDPM.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME HinderedSettlingDynamicsDPMFuncTest COMMAND $<TARGET_FILE:HinderedSettlingDynamicsDPM> --funcTest PROCESSES 4 LABELS longrun CONFIGURATIONS RelWithDbgInfo )

###################################################################################################
# Geometry tests
###################################################################################################

waLBerla_compile_test( FILES geometry/PeIntersectionRatioTest.cpp DEPENDS pe )
waLBerla_execute_test( NAME PeIntersectionRatioTest COMMAND $<TARGET_FILE:PeIntersectionRatioTest> PROCESSES 1 )

###################################################################################################
# Utility tests
###################################################################################################

waLBerla_compile_test( FILES utility/BodiesForceTorqueContainerTest.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME BodiesForceTorqueContainerTest COMMAND $<TARGET_FILE:BodiesForceTorqueContainerTest> PROCESSES 1 )
waLBerla_execute_test( NAME BodiesForceTorqueContainerParallelTest COMMAND $<TARGET_FILE:BodiesForceTorqueContainerTest> PROCESSES 3 )

waLBerla_compile_test( FILES utility/PeSubCyclingTest.cpp DEPENDS blockforest pe timeloop )
waLBerla_execute_test( NAME PeSubCyclingTest COMMAND $<TARGET_FILE:PeSubCyclingTest> PROCESSES 1 )
waLBerla_execute_test( NAME PeSubCyclingParallelTest COMMAND $<TARGET_FILE:PeSubCyclingTest> PROCESSES 3 )

