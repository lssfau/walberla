###################################################################################################
#
# Tests for field module
#
###################################################################################################

waLBerla_compile_test( FILES AccuracyEvaluationTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME AccuracyEvaluationTest4 COMMAND $<TARGET_FILE:AccuracyEvaluationTest> PROCESSES 4  )

waLBerla_compile_test( FILES AddToStorageTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME AddToStorageTest )

waLBerla_compile_test( FILES communication/FieldPackInfoTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME  FieldPackInfoTest )

waLBerla_compile_test( FILES distributors/DistributionTest.cpp)
waLBerla_execute_test( NAME DistributionTest )

waLBerla_compile_test( FILES FieldTest.cpp )
waLBerla_execute_test( NAME FieldTest )

waLBerla_compile_test( FILES FieldOfCustomTypesTest.cpp  )
waLBerla_execute_test( NAME FieldOfCustomTypesTest )

waLBerla_compile_test( FILES FieldTiming.cpp )
waLBerla_execute_test( NAME FieldTiming  )

waLBerla_compile_test( FILES FlagFieldTest.cpp)
waLBerla_execute_test( NAME FlagFieldTest )

waLBerla_compile_test( FILES StabilityCheckerTest.cpp DEPENDS blockforest field timeloop )
waLBerla_execute_test( NAME StabilityCheckerTest )

waLBerla_compile_test( FILES interpolators/InterpolationTest.cpp)
waLBerla_execute_test( NAME InterpolationTest )

waLBerla_compile_test( FILES interpolators/FieldInterpolationTest.cpp)
waLBerla_execute_test( NAME FieldInterpolationTest )

waLBerla_compile_test( FILES adaptors/AdaptorTest.cpp DEPENDS blockforest gui lbm )
waLBerla_execute_test( NAME AdaptorTest )

waLBerla_compile_test( FILES FieldGatherTest.cpp DEPENDS blockforest gui )
waLBerla_execute_test( NAME FieldGatherTest3 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 3 DEPENDS_ON_TARGETS FieldGatherTest )
waLBerla_execute_test( NAME FieldGatherTest1 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 1 DEPENDS_ON_TARGETS FieldGatherTest )

waLBerla_compile_test( FILES FieldFileIOTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME FieldFileIOTest1Proc  COMMAND $<TARGET_FILE:FieldFileIOTest> PROCESSES 1  )
waLBerla_execute_test( NAME FieldFileIOTest2Proc  COMMAND $<TARGET_FILE:FieldFileIOTest> PROCESSES 2  )
waLBerla_execute_test( NAME FieldFileIOTest4Proc  COMMAND $<TARGET_FILE:FieldFileIOTest> PROCESSES 4  )
waLBerla_execute_test( NAME FieldFileIOTest16Proc COMMAND $<TARGET_FILE:FieldFileIOTest> PROCESSES 16 )
#serialize runs of tests to avoid i/o conflicts when running ctest with -jN
if( WALBERLA_BUILD_WITH_MPI )
   set_property( TEST FieldFileIOTest2Proc  PROPERTY DEPENDS FieldFileIOTest1Proc )
   set_property( TEST FieldFileIOTest4Proc  PROPERTY DEPENDS FieldFileIOTest2Proc )
   set_property( TEST FieldFileIOTest16Proc PROPERTY DEPENDS FieldFileIOTest4Proc )
endif( WALBERLA_BUILD_WITH_MPI )

if( WALBERLA_BUILD_WITH_MPI )
   waLBerla_compile_test( FILES FieldMPIDatatypesTest.cpp DEPENDS blockforest )
   waLBerla_execute_test( NAME FieldMPIDatatypesTestRelease  COMMAND $<TARGET_FILE:FieldMPIDatatypesTest> PROCESSES 1                CONFIGURATIONS Release RelWithDbgInfo )
   waLBerla_execute_test( NAME FieldMPIDatatypesTestDebug    COMMAND $<TARGET_FILE:FieldMPIDatatypesTest> PROCESSES 1 LABELS longrun CONFIGURATIONS Debug DebugOptimized   )
endif( WALBERLA_BUILD_WITH_MPI )



# CodeGen Tests
if( WALBERLA_BUILD_WITH_CODEGEN )
waLBerla_generate_target_from_python(NAME CodegenJacobiCPUGeneratedJacobiKernel FILE codegen/JacobiKernel.py
      OUT_FILES JacobiKernel2D.cpp JacobiKernel2D.h JacobiKernel3D.cpp JacobiKernel3D.h )
waLBerla_compile_test( FILES codegen/CodegenJacobiCPU.cpp DEPENDS gui timeloop CodegenJacobiCPUGeneratedJacobiKernel)
waLBerla_execute_test( NAME CodegenJacobiCPU )

waLBerla_generate_target_from_python(NAME SweepCollectionKernel FILE codegen/SweepCollection.py
        OUT_FILES SweepCollection.h SweepCollection.cpp)
waLBerla_compile_test( FILES codegen/SweepCollection.cpp DEPENDS timeloop SweepCollectionKernel)
waLBerla_execute_test( NAME SweepCollection )

waLBerla_generate_target_from_python(NAME CodegenPoissonCPUGeneratedKernel FILE codegen/Poisson.py
      OUT_FILES Poisson.cpp Poisson.h )
waLBerla_compile_test( FILES codegen/CodegenPoissonCPU.cpp DEPENDS gui timeloop CodegenPoissonCPUGeneratedKernel)
waLBerla_execute_test( NAME CodegenPoissonCPU )

waLBerla_generate_target_from_python(NAME CodeGenMultipleFieldSwaps FILE codegen/MultipleFieldSwaps.py
        OUT_FILES MultipleFieldSwaps.cpp MultipleFieldSwaps.h )
waLBerla_compile_test( FILES codegen/MultipleFieldSwaps.cpp DEPENDS gui timeloop CodeGenMultipleFieldSwaps)
waLBerla_execute_test( NAME MultipleFieldSwaps )

waLBerla_generate_target_from_python(NAME CodegenGeneratedCPUFieldPackInfo FILE codegen/GeneratedFieldPackInfoTest.py
        OUT_FILES ScalarFieldCommunication.cpp ScalarFieldCommunication.h
                  ScalarFieldPullReduction.cpp ScalarFieldPullReduction.h )
waLBerla_compile_test( FILES codegen/GeneratedFieldPackInfoTest.cpp
        DEPENDS blockforest core field CodegenGeneratedCPUFieldPackInfo )
waLBerla_execute_test( NAME GeneratedFieldPackInfoTest )

waLBerla_generate_target_from_python(NAME CodeGenEK FILE codegen/EK.py
        OUT_FILES EKFlux.cpp EKFlux.h EKContinuity.cpp EKContinuity.h )
waLBerla_compile_test( FILES codegen/EK.cpp DEPENDS blockforest core field timeloop CodeGenEK)
waLBerla_execute_test( NAME EK )
endif()
