#############################################################################################################################
#
# Tests for walberla::lbm module
#
#############################################################################################################################

waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_test_executable ( SweepEquivalenceTest SweepEquivalenceTest.cpp )
target_link_libraries ( SweepEquivalenceTest walberla::blockforest walberla::timeloop walberla::core walberla::field
        walberla::lbm)
waLBerla_execute_test( NAME SweepEquivalenceTest )

waLBerla_add_test_executable ( BoundaryHandlingCommunication BoundaryHandlingCommunication.cpp )
target_link_libraries ( BoundaryHandlingCommunication walberla::blockforest walberla::timeloop walberla::core
        walberla::lbm )
waLBerla_execute_test( NAME BoundaryHandlingCommunication PROCESSES 8 )

waLBerla_add_test_executable ( UnrollTest UnrollTest.cpp )
target_link_libraries ( UnrollTest walberla::core walberla::core walberla::lbm )


waLBerla_add_test_executable ( SimplePABTest boundary/SimplePABTest.cpp )
target_link_libraries ( SimplePABTest walberla::field walberla::blockforest walberla::timeloop walberla::vtk
        walberla::core walberla::lbm )


waLBerla_add_test_executable ( SimpleDiffusionDirichlet boundary/SimpleDiffusionDirichlet.cpp )
target_link_libraries ( SimpleDiffusionDirichlet walberla::field walberla::blockforest walberla::geometry
        walberla::timeloop walberla::vtk walberla::core walberla::lbm )
waLBerla_execute_test( NAME SimpleDiffusionDirichletDef COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 4 -w 2 -o 1.2 -d 3 -t 333 -c 0 -r 0 )
waLBerla_execute_test( NAME SimpleDiffusionDirichletRef COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 4 -w 4 -o 1.2 -d 3 -t 3   -c 0 -r 1 )
waLBerla_execute_test( NAME SimpleDiffusionDirichletCav COMMAND $<TARGET_FILE:SimpleDiffusionDirichlet>  -l 5 -w 3 -o 1.2 -d 3 -t 333 -c 1 )

waLBerla_add_test_executable ( BoundaryForce boundary/BoundaryForce.cpp )
target_link_libraries ( BoundaryForce walberla::blockforest walberla::timeloop walberla::core walberla::lbm )
waLBerla_execute_test( NAME BoundaryForceShort
        COMMAND $<TARGET_FILE:BoundaryForce> 10 PROCESSES 2 CONFIGURATIONS Debug DebugOptimized )
waLBerla_execute_test( NAME BoundaryForce
        COMMAND $<TARGET_FILE:BoundaryForce> PROCESSES 2 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( BoundaryForceCouette boundary/BoundaryForceCouette.cpp )
target_link_libraries ( BoundaryForceCouette walberla::blockforest walberla::timeloop walberla::core walberla::lbm )
waLBerla_execute_test( NAME BoundaryForceCouetteShort
        COMMAND $<TARGET_FILE:BoundaryForceCouette> 10 PROCESSES 2 CONFIGURATIONS Debug DebugOptimized )
waLBerla_execute_test( NAME BoundaryForceCouette
        COMMAND $<TARGET_FILE:BoundaryForceCouette> PROCESSES 2 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( DiffusionDirichlet boundary/DiffusionDirichlet.cpp )
target_link_libraries ( DiffusionDirichlet walberla::field walberla::blockforest walberla::geometry walberla::timeloop
        walberla::vtk walberla::core walberla::lbm )
waLBerla_execute_test( NAME DiffusionDirichletTest1
        COMMAND $<TARGET_FILE:DiffusionDirichlet> -l 16 -w 1 -o 1.2  -v 0.1  -e 0.097711 -t  50 )
waLBerla_execute_test( NAME DiffusionDirichletTest2
        COMMAND $<TARGET_FILE:DiffusionDirichlet> -l 16 -w 1 -o 1.79 -v 0.05 -e 0.295170 -t 100 )


waLBerla_add_test_executable ( DiffusionTest DiffusionTest.cpp )
target_link_libraries ( DiffusionTest walberla::field walberla::blockforest walberla::timeloop walberla::vtk
        walberla::postprocessing walberla::core walberla::lbm )
waLBerla_execute_test( NAME  DiffusionTestCorr
        COMMAND $<TARGET_FILE:DiffusionTest> -c 1  -dx 0.01  -dt 0.002  -d 0.005  -v 5  -err 0.0025
        DEPENDS_ON_TARGETS DiffusionTest )
waLBerla_execute_test( NAME  DiffusionTestNone
        COMMAND $<TARGET_FILE:DiffusionTest> -c 0  -dx 0.01  -dt 0.002  -d 0.005  -v 5  -err 0.0040
        DEPENDS_ON_TARGETS DiffusionTest )


waLBerla_add_test_executable ( Uniformity refinement/Uniformity.cpp )
target_link_libraries ( Uniformity walberla::blockforest walberla::stencil walberla::core walberla::lbm )
waLBerla_execute_test( NAME UniformityShortTest COMMAND $<TARGET_FILE:Uniformity> --shortrun PROCESSES 4 )
waLBerla_execute_test( NAME UniformityLongTest
        COMMAND $<TARGET_FILE:Uniformity> PROCESSES 4 LABELS longrun verylongrun CONFIGURATIONS Release RelWithDbgInfo )


waLBerla_add_test_executable ( CommunicationEquivalence refinement/CommunicationEquivalence.cpp )
target_link_libraries ( CommunicationEquivalence walberla::blockforest walberla::stencil walberla::core walberla::lbm )
waLBerla_execute_test( NAME CommunicationEquivalenceShortTest
        COMMAND $<TARGET_FILE:CommunicationEquivalence> --shortrun PROCESSES 4 )
waLBerla_execute_test( NAME CommunicationEquivalenceLongTest
        COMMAND $<TARGET_FILE:CommunicationEquivalence> PROCESSES 4 LABELS longrun
        CONFIGURATIONS Release RelWithDbgInfo )


waLBerla_add_test_executable ( NonConstantDiffusion refinement/NonConstantDiffusion.cpp )
target_link_libraries ( NonConstantDiffusion walberla::field walberla::blockforest walberla::geometry walberla::timeloop
        walberla::vtk walberla::core walberla::lbm )


waLBerla_add_test_executable ( IntersectionRatioTest geometry/IntersectionRatioTest.cpp )
target_link_libraries ( IntersectionRatioTest walberla::geometry walberla::core walberla::lbm )
waLBerla_execute_test( NAME IntersectionRatioTest COMMAND $<TARGET_FILE:IntersectionRatioTest> )

waLBerla_add_test_executable ( PermeabilityTest evaluations/PermeabilityTest.cpp )
target_link_libraries ( PermeabilityTest walberla::field walberla::blockforest walberla::geometry walberla::timeloop
        walberla::vtk walberla::core walberla::lbm )
waLBerla_execute_test( NAME PermeabilityTest_TRT_16_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 16 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.005 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_32_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 32 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.005 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_1 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 1 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_4 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )
waLBerla_execute_test( NAME PermeabilityTest_TRT_64_8 COMMAND $<TARGET_FILE:PermeabilityTest> --length 64 --kappa 0.85 --omega 0.3 --collisionModel TRT --timesteps 1001 --epsilon 0.025 PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( UniformPdfFieldInitializerTest initializer/UniformPdfFieldInitializerTest.cpp )
target_link_libraries( UniformPdfFieldInitializerTest walberla::lbm )
waLBerla_execute_test( NAME UniformPdfFieldInitializerTest COMMAND $<TARGET_FILE:UniformPdfFieldInitializerTest> ${CMAKE_CURRENT_SOURCE_DIR}/PdfFieldInitializerTest.prm PROCESSES 4 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( NonUniformPdfFieldInitializerTest initializer/NonUniformPdfFieldInitializerTest.cpp )
target_link_libraries( NonUniformPdfFieldInitializerTest walberla::lbm )
waLBerla_execute_test( NAME NonUniformPdfFieldInitializerTest CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( QCriterionTest field/QCriterionTest.cpp )
target_link_libraries ( QCriterionTest walberla::field walberla::blockforest walberla::core )
waLBerla_execute_test( NAME QCriterionTest )

add_subdirectory(field)
waLBerla_add_test_executable ( AlignedPDFfieldTest field/AlignedPDFfieldTest.cpp )
target_link_libraries ( AlignedPDFfieldTest walberla::field walberla::blockforest walberla::lbm walberla::core )
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
            GeneratedOutflowBC.h)
    waLBerla_add_test_executable ( GeneratedOutflowBC codegen/GeneratedOutflowBC.cpp )
    target_link_libraries ( GeneratedOutflowBC GeneratedOutflowBCGenerated walberla::core
            walberla::blockforest walberla::field walberla::geometry walberla::timeloop )
    waLBerla_execute_test( NAME GeneratedOutflowBC COMMAND $<TARGET_FILE:GeneratedOutflowBC> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/GeneratedOutflowBC.prm  )

    waLBerla_generate_target_from_python(NAME GeneratedFreeSlipGenerated
            FILE codegen/GeneratedFreeSlip.py
            OUT_FILES GeneratedFreeSlip_Sweep.cpp GeneratedFreeSlip_Sweep.h
            GeneratedFreeSlip_MacroSetter.cpp GeneratedFreeSlip_MacroSetter.h
            GeneratedFreeSlip_NoSlip.cpp GeneratedFreeSlip_NoSlip.h
            GeneratedFreeSlip_FreeSlip.cpp GeneratedFreeSlip_FreeSlip.h
            GeneratedFreeSlip_PackInfoEven.cpp GeneratedFreeSlip_PackInfoEven.h
            GeneratedFreeSlip_PackInfoOdd.cpp GeneratedFreeSlip_PackInfoOdd.h
            GeneratedFreeSlip.h)
    waLBerla_add_test_executable ( GeneratedFreeSlip codegen/GeneratedFreeSlip.cpp )
    target_link_libraries ( GeneratedFreeSlip GeneratedFreeSlipGenerated walberla::core
            walberla::blockforest walberla::field walberla::geometry walberla::timeloop )
    waLBerla_execute_test( NAME GeneratedFreeSlip COMMAND $<TARGET_FILE:GeneratedFreeSlip> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/GeneratedFreeSlip.prm  )


    waLBerla_generate_target_from_python(NAME LbCodeGenerationExampleGenerated
            FILE codegen/LbCodeGenerationExample.py
            OUT_FILES LbCodeGenerationExample_LatticeModel.cpp LbCodeGenerationExample_LatticeModel.h
            LbCodeGenerationExample_NoSlip.cpp LbCodeGenerationExample_NoSlip.h
            LbCodeGenerationExample_UBB.cpp LbCodeGenerationExample_UBB.h
            LbCodeGenerationExample.h)
    waLBerla_add_test_executable ( LbCodeGenerationExample codegen/LbCodeGenerationExample.cpp )
    target_link_libraries ( LbCodeGenerationExample LbCodeGenerationExampleGenerated walberla::core walberla::lbm )

    waLBerla_generate_target_from_python(NAME FluctuatingMRTGenerated FILE codegen/FluctuatingMRT.py
            OUT_FILES FluctuatingMRT_LatticeModel.cpp FluctuatingMRT_LatticeModel.h )
    waLBerla_add_test_executable ( FluctuatingMRT codegen/FluctuatingMRT.cpp )
    target_link_libraries ( FluctuatingMRT FluctuatingMRTGenerated walberla::core walberla::lbm )
    waLBerla_execute_test( NAME FluctuatingMRT COMMAND $<TARGET_FILE:FluctuatingMRT> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/FluctuatingMRT.prm PROCESSES 1)

    waLBerla_generate_target_from_python(NAME FieldLayoutAndVectorizationTestGenerated FILE codegen/FieldLayoutAndVectorizationTest.py
            OUT_FILES FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel.cpp FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel.h
            FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel.cpp FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel.h
            FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel.cpp FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel.h
            FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel.cpp FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel.h )

    waLBerla_add_test_executable ( FieldLayoutAndVectorizationTest codegen/FieldLayoutAndVectorizationTest.cpp )
    target_link_libraries ( FieldLayoutAndVectorizationTest FieldLayoutAndVectorizationTestGenerated walberla::core
            walberla::lbm walberla::blockforest )
    waLBerla_execute_test( NAME FieldLayoutAndVectorizationTest )

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

    waLBerla_add_test_executable ( InplaceStreamingCodegenTest codegen/InplaceStreamingCodegenTest.cpp )
    target_link_libraries ( InplaceStreamingCodegenTest InplaceStreamingCodegen walberla::core walberla::lbm )

    waLBerla_execute_test( NAME InplaceStreamingCodegenTest
            COMMAND $<TARGET_FILE:InplaceStreamingCodegenTest> ${CMAKE_CURRENT_SOURCE_DIR}/codegen/InplaceStreamingCodegen3D.prm  )


    waLBerla_generate_target_from_python( NAME StreamInCellIntervalCodegen
            FILE codegen/StreamInCellIntervalCodegen.py
            OUT_FILES GeneratedLatticeModel.h GeneratedLatticeModel.cpp)

    waLBerla_add_test_executable ( StreamInCellIntervalCodegenTest codegen/StreamInCellIntervalCodegenTest.cpp )
    target_link_libraries ( StreamInCellIntervalCodegenTest walberla::blockforest walberla::field walberla::lbm
            walberla::timeloop StreamInCellIntervalCodegen walberla::core )

    waLBerla_execute_test( NAME     StreamInCellIntervalCodegenTest )

endif()

# Free Surface
waLBerla_add_test_executable ( BubbleInitializationTest free_surface/bubble_model/BubbleInitializationTest.cpp )
target_link_libraries ( BubbleInitializationTest walberla::blockforest walberla::field walberla::geometry walberla::lbm
        walberla::timeloop walberla::core )
waLBerla_execute_test ( NAME BubbleInitializationTest PROCESSES 2 )

file ( COPY free_surface/bubble_model/MergeAndSplitTestUnconnected.png
        free_surface/bubble_model/MergeAndSplitTestConnected.png
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
waLBerla_add_test_executable ( MergeAndSplitTest
        free_surface/bubble_model/MergeAndSplitTest.cpp
        free_surface/bubble_model/BubbleBodyMover.impl.h
        free_surface/bubble_model/BubbleModelTester.impl.h
)
target_link_libraries ( MergeAndSplitTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop )
waLBerla_execute_test ( NAME MergeAndSplitTest PROCESSES 10 )

waLBerla_add_test_executable ( MergeInformationTest free_surface/bubble_model/MergeInformationTest.cpp )
target_link_libraries ( MergeInformationTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop walberla::core )
waLBerla_execute_test( NAME      MergeInformationTest
        PROCESSES 3 )

waLBerla_add_test_executable ( MovingSpheresTest free_surface/bubble_model/MovingSpheresTest.cpp
        free_surface/bubble_model/BubbleBodyMover.impl.h
        free_surface/bubble_model/BubbleModelTester.impl.h
)
target_link_libraries ( MovingSpheresTest walberla::blockforest walberla::field walberla::geometry walberla::lbm
        walberla::timeloop )

waLBerla_add_test_executable ( RegionalFloodFillTest free_surface/bubble_model/RegionalFloodFillTest.cpp )
target_link_libraries ( RegionalFloodFillTest walberla::field walberla::lbm walberla::core )
waLBerla_execute_test( NAME     RegionalFloodFillTest )

waLBerla_add_test_executable ( SplitDetectionTest free_surface/bubble_model/SplitDetectionTest.cpp )
target_link_libraries ( SplitDetectionTest walberla::field walberla::lbm walberla::core )
waLBerla_execute_test( NAME     SplitDetectionTest )

waLBerla_add_test_executable ( AdvectionTest free_surface/dynamics/AdvectionTest.cpp )
target_link_libraries ( AdvectionTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop walberla::core )
waLBerla_execute_test( NAME     AdvectionTest )

waLBerla_add_test_executable ( CellConversionTest free_surface/dynamics/CellConversionTest.cpp )
target_link_libraries ( CellConversionTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop walberla::core )
waLBerla_execute_test( NAME     CellConversionTest )

if( WALBERLA_BUILD_WITH_CODEGEN )
    walberla_generate_target_from_python( NAME      LatticeModelGenerationFreeSurfacePython
            FILE      free_surface/dynamics/LatticeModelGenerationFreeSurface.py
            OUT_FILES GeneratedLatticeModel_FreeSurface.cpp
            GeneratedLatticeModel_FreeSurface.h )
    waLBerla_add_test_executable ( CodegenTest free_surface/dynamics/CodegenTest.cpp )
    target_link_libraries ( CodegenTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
            LatticeModelGenerationFreeSurfacePython walberla::core )
    waLBerla_execute_test( NAME     CodegenTest )
endif()

waLBerla_add_test_executable ( ExcessMassDistributionFallbackTest free_surface/dynamics/ExcessMassDistributionFallbackTest.cpp )
target_link_libraries ( ExcessMassDistributionFallbackTest walberla::blockforest walberla::field walberla::lbm
        walberla::timeloop walberla::core )
waLBerla_execute_test( NAME     ExcessMassDistributionFallbackTest )

waLBerla_add_test_executable ( ExcessMassDistributionParallelTest free_surface/dynamics/ExcessMassDistributionParallelTest.cpp )
target_link_libraries ( ExcessMassDistributionParallelTest walberla::blockforest walberla::field walberla::lbm
        walberla::timeloop walberla::core )
waLBerla_execute_test( NAME ExcessMassDistributionParallelTest PROCESSES 2)

waLBerla_add_test_executable ( InflowTest free_surface/dynamics/InflowTest.cpp )
target_link_libraries ( InflowTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     InflowTest )

waLBerla_add_test_executable ( LoadBalancingTest free_surface/LoadBalancingTest.cpp )
target_link_libraries ( LoadBalancingTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME         LoadBalancingTest
        PROCESSES    4)

waLBerla_add_test_executable ( PdfReconstructionFreeSlipTest free_surface/dynamics/PdfReconstructionFreeSlipTest.cpp )
target_link_libraries ( PdfReconstructionFreeSlipTest walberla::blockforest walberla::field walberla::lbm
        walberla::timeloop walberla::core )
waLBerla_execute_test( NAME     PdfReconstructionFreeSlipTest )

waLBerla_add_test_executable ( PdfReconstructionTest free_surface/dynamics/PdfReconstructionTest.cpp )
target_link_libraries ( PdfReconstructionTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     PdfReconstructionTest )

waLBerla_add_test_executable ( PdfRefillingTest free_surface/dynamics/PdfRefillingTest.cpp )
target_link_libraries ( PdfRefillingTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     PdfRefillingTest )

waLBerla_add_test_executable ( CellFluidVolumeTest free_surface/surface_geometry/CellFluidVolumeTest.cpp )
target_link_libraries ( CellFluidVolumeTest walberla::lbm walberla::core )
waLBerla_execute_test( NAME     CellFluidVolumeTest )

waLBerla_add_test_executable ( CurvatureOfSineTest free_surface/surface_geometry/CurvatureOfSineTest.cpp )
target_link_libraries ( CurvatureOfSineTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )

waLBerla_add_test_executable ( WettingConversionTest free_surface/dynamics/WettingConversionTest.cpp )
target_link_libraries ( WettingConversionTest walberla::blockforest walberla::field walberla::geometry walberla::lbm
        walberla::timeloop walberla::core )

waLBerla_add_test_executable ( CurvatureOfSphereTest free_surface/surface_geometry/CurvatureOfSphereTest.cpp )
target_link_libraries ( CurvatureOfSphereTest walberla::blockforest walberla::field walberla::geometry walberla::lbm
        walberla::timeloop walberla::core )

waLBerla_add_test_executable ( DetectWettingTest free_surface/surface_geometry/DetectWettingTest.cpp )
target_link_libraries ( DetectWettingTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     DetectWettingTest )

waLBerla_add_test_executable ( GetInterfacePointTest free_surface/surface_geometry/GetInterfacePointTest.cpp )
target_link_libraries ( GetInterfacePointTest walberla::lbm walberla::core )
waLBerla_execute_test( NAME     GetInterfacePointTest )

waLBerla_add_test_executable ( NormalsOfSineTest free_surface/surface_geometry/NormalsOfSineTest.cpp )
target_link_libraries ( NormalsOfSineTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )

waLBerla_add_test_executable ( NormalsOfSphereTest free_surface/surface_geometry/NormalsOfSphereTest.cpp )
target_link_libraries ( NormalsOfSphereTest walberla::blockforest walberla::field walberla::geometry walberla::lbm
        walberla::timeloop walberla::core )

waLBerla_add_test_executable ( NormalsEquivalenceTest free_surface/surface_geometry/NormalsEquivalenceTest.cpp )
target_link_libraries ( NormalsEquivalenceTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     NormalsEquivalenceTest )

waLBerla_add_test_executable ( NormalsNearSolidTest free_surface/surface_geometry/NormalsNearSolidTest.cpp )
target_link_libraries ( NormalsNearSolidTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     NormalsNearSolidTest )

waLBerla_add_test_executable ( ObstacleFillLevelTest free_surface/surface_geometry/ObstacleFillLevelTest.cpp )
target_link_libraries ( ObstacleFillLevelTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     ObstacleFillLevelTest )

waLBerla_add_test_executable ( ObstacleNormalsTest free_surface/surface_geometry/ObstacleNormalsTest.cpp )
target_link_libraries ( ObstacleNormalsTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     ObstacleNormalsTest )

waLBerla_add_test_executable ( WettingCurvatureTest free_surface/surface_geometry/WettingCurvatureTest.cpp )
target_link_libraries ( WettingCurvatureTest walberla::blockforest walberla::field walberla::lbm walberla::timeloop
        walberla::core )
waLBerla_execute_test( NAME     WettingCurvatureTest )

#     Slow-running tests

waLBerla_execute_test( NAME     MovingSpheresTest PROCESSES 2 LABELS longrun )
waLBerla_execute_test( NAME     WettingConversionTest LABELS longrun )
waLBerla_execute_test( NAME     CurvatureOfSineTest LABELS longrun )
waLBerla_execute_test( NAME     CurvatureOfSphereTest LABELS longrun )
waLBerla_execute_test( NAME     NormalsOfSineTest LABELS longrun )
waLBerla_execute_test( NAME     NormalsOfSphereTest LABELS longrun )

