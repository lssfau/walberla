#############################################################################################################################
#
# Tests for lbm module
#
#############################################################################################################################

waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_compile_test( FILES SweepEquivalenceTest.cpp DEPENDS blockforest timeloop )
waLBerla_execute_test( NAME SweepEquivalenceTest )

waLBerla_compile_test( FILES BoundaryHandlingCommunication.cpp DEPENDS blockforest timeloop )
waLBerla_execute_test( NAME BoundaryHandlingCommunication PROCESSES 8 )

waLBerla_compile_test( FILES UnrollTest.cpp  )


waLBerla_compile_test( FILES boundary/SimplePABTest.cpp DEPENDS field blockforest timeloop vtk )


waLBerla_compile_test( FILES boundary/SimpleDiffusionDirichlet.cpp DEPENDS field blockforest geometry timeloop vtk )
waLBerla_execute_test( NAME SimpleDiffusionDirichletDef COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 4 -w 2 -o 1.2 -d 3 -t 333 -c 0 -r 0 )
waLBerla_execute_test( NAME SimpleDiffusionDirichletRef COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 4 -w 4 -o 1.2 -d 3 -t 3   -c 0 -r 1 )
waLBerla_execute_test( NAME SimpleDiffusionDirichletCav COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 5 -w 3 -o 1.2 -d 3 -t 333 -c 1 )

waLBerla_compile_test( FILES boundary/BoundaryForce.cpp DEPENDS blockforest timeloop )
waLBerla_execute_test( NAME BoundaryForceShort COMMAND $<TARGET_FILE:BoundaryForce> 10 PROCESSES 2 CONFIGURATIONS Debug DebugOptimized )
waLBerla_execute_test( NAME BoundaryForce COMMAND COMMAND $<TARGET_FILE:BoundaryForce> PROCESSES 2 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES boundary/BoundaryForceCouette.cpp DEPENDS blockforest timeloop )
waLBerla_execute_test( NAME BoundaryForceCouetteShort COMMAND $<TARGET_FILE:BoundaryForceCouette> 10 PROCESSES 2 CONFIGURATIONS Debug DebugOptimized )
waLBerla_execute_test( NAME BoundaryForceCouette COMMAND COMMAND $<TARGET_FILE:BoundaryForceCouette> PROCESSES 2 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES boundary/DiffusionDirichlet.cpp DEPENDS field blockforest geometry timeloop vtk )
waLBerla_execute_test( NAME DiffusionDirichletTest1 COMMAND $<TARGET_FILE:DiffusionDirichlet>  -l 16 -w 1 -o 1.2  -v 0.1  -e 0.097711 -t  50 )
waLBerla_execute_test( NAME DiffusionDirichletTest2 COMMAND $<TARGET_FILE:DiffusionDirichlet>  -l 16 -w 1 -o 1.79 -v 0.05 -e 0.295170 -t 100 )


waLBerla_compile_test( FILES DiffusionTest.cpp DEPENDS field blockforest timeloop vtk gui postprocessing)
waLBerla_execute_test( NAME  DiffusionTestCorr COMMAND $<TARGET_FILE:DiffusionTest>  -c 1  -dx 0.01  -dt 0.002  -d 0.005  -v 5  -err 0.0025 )
waLBerla_execute_test( NAME  DiffusionTestNone COMMAND $<TARGET_FILE:DiffusionTest>  -c 0  -dx 0.01  -dt 0.002  -d 0.005  -v 5  -err 0.0040 )


waLBerla_compile_test( FILES refinement/Uniformity.cpp DEPENDS blockforest stencil )
waLBerla_execute_test( NAME UniformityShortTest COMMAND $<TARGET_FILE:Uniformity> --shortrun PROCESSES 4                )
waLBerla_execute_test( NAME UniformityLongTest  COMMAND $<TARGET_FILE:Uniformity>            PROCESSES 4 LABELS longrun verylongrun CONFIGURATIONS Release RelWithDbgInfo )


waLBerla_compile_test( FILES refinement/CommunicationEquivalence.cpp DEPENDS blockforest stencil )
waLBerla_execute_test( NAME CommunicationEquivalenceShortTest COMMAND $<TARGET_FILE:CommunicationEquivalence> --shortrun PROCESSES 4                )
waLBerla_execute_test( NAME CommunicationEquivalenceLongTest  COMMAND $<TARGET_FILE:CommunicationEquivalence>            PROCESSES 4 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )


waLBerla_compile_test( FILES refinement/NonConstantDiffusion.cpp DEPENDS field blockforest geometry timeloop vtk )


waLBerla_compile_test( FILES geometry/IntersectionRatioTest.cpp DEPENDS geometry )
waLBerla_execute_test( NAME IntersectionRatioTest COMMAND $<TARGET_FILE:IntersectionRatioTest> )

waLBerla_compile_test( FILES evaluations/PermeabilityTest.cpp DEPENDS field blockforest geometry timeloop vtk )
waLBerla_execute_test( NAME PermeabilityTest_TRT_16_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 16 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.005 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_32_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 32 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.005 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_1 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 1 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_8 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES initializer/PdfFieldInitializerTest.cpp )
waLBerla_execute_test( NAME PdfFieldInitializerTest COMMAND $<TARGET_FILE:PdfFieldInitializerTest> ${CMAKE_CURRENT_SOURCE_DIR}/PdfFieldInitializerTest.prm PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES SuViscoelasticityTest.cpp DEPENDS field blockforest timeloop)
waLBerla_execute_test( NAME  SuViscoelasticityLongTest COMMAND $<TARGET_FILE:SuViscoelasticityTest> ${CMAKE_CURRENT_SOURCE_DIR}/Su.prm LABELS longrun)
waLBerla_execute_test( NAME  SuViscoelasticityShortTest COMMAND $<TARGET_FILE:SuViscoelasticityTest> ${CMAKE_CURRENT_SOURCE_DIR}/Su.prm -Parameters.shortrun=true )

waLBerla_compile_test( FILES field/QCriterionTest.cpp DEPENDS field blockforest )
waLBerla_execute_test( NAME QCriterionTest )

add_subdirectory(field)
waLBerla_compile_test( FILES field/AlignedPDFfieldTest.cpp DEPENDS field blockforest lbm )
waLBerla_execute_test( NAME AlignedPDFfieldTest COMMAND $<TARGET_FILE:AlignedPDFfieldTest> ${CMAKE_CURRENT_SOURCE_DIR}/field/AlignedPDFfieldTest.prm )

# Code Generation
if( WALBERLA_BUILD_WITH_CODEGEN )
add_subdirectory(codegen)

waLBerla_generate_target_from_python(NAME GeneratedOutflowBCGenerated
        FILE codegen/GeneratedOutflowBC.py
        OUT_FILES GeneratedOutflowBC_Sweep.cpp GeneratedOutflowBC_Sweep.h
        GeneratedOutflowBC_MacroSetter.cpp GeneratedOutflowBC_MacroSetter.h
        GeneratedOutflowBC_Dynamic_UBB.cpp GeneratedOutflowBC_Dynamic_UBB.h
        GeneratedOutflowBC_Static_UBB.cpp GeneratedOutflowBC_Static_UBB.h
        GeneratedOutflowBC_NoSlip.cpp GeneratedOutflowBC_NoSlip.h
        GeneratedOutflowBC_Outflow.cpp GeneratedOutflowBC_Outflow.h
        GeneratedOutflowBC_PackInfo.cpp GeneratedOutflowBC_PackInfo.h
        GeneratedOutflowBC_InfoHeader.h)
waLBerla_compile_test( FILES codegen/GeneratedOutflowBC.cpp DEPENDS GeneratedOutflowBCGenerated)
waLBerla_execute_test( NAME GeneratedOutflowBC COMMAND $<TARGET_FILE:GeneratedOutflowBC> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/GeneratedOutflowBC.prm  )


waLBerla_generate_target_from_python(NAME LbCodeGenerationExampleGenerated
      FILE codegen/LbCodeGenerationExample.py
      OUT_FILES LbCodeGenerationExample_LatticeModel.cpp LbCodeGenerationExample_LatticeModel.h
      LbCodeGenerationExample_NoSlip.cpp LbCodeGenerationExample_NoSlip.h
      LbCodeGenerationExample_UBB.cpp LbCodeGenerationExample_UBB.h )
waLBerla_compile_test( FILES codegen/LbCodeGenerationExample.cpp DEPENDS LbCodeGenerationExampleGenerated)

waLBerla_generate_target_from_python(NAME FluctuatingMRTGenerated FILE codegen/FluctuatingMRT.py
                              OUT_FILES FluctuatingMRT_LatticeModel.cpp FluctuatingMRT_LatticeModel.h )
waLBerla_compile_test( FILES codegen/FluctuatingMRT.cpp DEPENDS FluctuatingMRTGenerated)
waLBerla_execute_test( NAME FluctuatingMRT COMMAND $<TARGET_FILE:FluctuatingMRT> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/FluctuatingMRT.prm PROCESSES 1)

waLBerla_generate_target_from_python(NAME FieldLayoutAndVectorizationTestGenerated FILE codegen/FieldLayoutAndVectorizationTest.py
                                     OUT_FILES FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel.cpp FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel.h
                                               FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel.cpp FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel.h
                                               FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel.cpp FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel.h
                                               FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel.cpp FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel.h )
waLBerla_compile_test( FILES codegen/FieldLayoutAndVectorizationTest.cpp DEPENDS FieldLayoutAndVectorizationTestGenerated)

waLBerla_generate_target_from_python(NAME LbmPackInfoGenerationTestCodegen FILE codegen/LbmPackInfoGenerationTest.py
                                     OUT_FILES AccessorBasedPackInfoEven.cpp AccessorBasedPackInfoEven.h
                                               AccessorBasedPackInfoOdd.cpp AccessorBasedPackInfoOdd.h
                                               FromKernelPackInfoPull.cpp FromKernelPackInfoPull.h
                                               FromKernelPackInfoPush.cpp FromKernelPackInfoPush.h)

waLBerla_link_files_to_builddir( "diff_packinfos.sh" )
waLBerla_execute_test( NAME LbmPackInfoGenerationDiffTest COMMAND bash diff_packinfos.sh )


waLBerla_generate_target_from_python( NAME InplaceStreamingCodegen
      FILE codegen/InplaceStreamingCodegen.py
      OUT_FILES PullSweep.h PullSweep.cpp
                PullNoSlip.h PullNoSlip.cpp
                PullUBB.h PullUBB.cpp
                PullOutflow.h PullOutflow.cpp
                PullInit.h PullInit.cpp
                InPlaceSweep.h InPlaceSweep.cpp
                InPlaceNoSlip.h InPlaceNoSlip.cpp
                InPlaceUBB.h InPlaceUBB.cpp
                InPlaceOutflow.h InPlaceOutflow.cpp
                InPlaceInit.h InPlaceInit.cpp
                InplaceStreamingCodegen.h)

waLBerla_compile_test( FILES codegen/InplaceStreamingCodegenTest.cpp
                       DEPENDS InplaceStreamingCodegen )

waLBerla_execute_test( NAME InplaceStreamingCodegenTest
                       COMMAND $<TARGET_FILE:InplaceStreamingCodegenTest> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/InplaceStreamingCodegen3D.prm  )

endif()
