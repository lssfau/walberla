###################################################################################################
#
# Tests for walberla::field module
#
###################################################################################################

waLBerla_add_test_executable ( AccuracyEvaluationTest AccuracyEvaluationTest.cpp )
target_link_libraries ( AccuracyEvaluationTest walberla::blockforest walberla::core walberla::field )
waLBerla_execute_test( NAME AccuracyEvaluationTest4 COMMAND $<TARGET_FILE:AccuracyEvaluationTest> PROCESSES 4  )

waLBerla_add_test_executable ( AddToStorageTest AddToStorageTest.cpp )
target_link_libraries ( AddToStorageTest walberla::blockforest walberla::core walberla::field )
waLBerla_execute_test( NAME AddToStorageTest )

waLBerla_add_test_executable ( FieldPackInfoTest communication/FieldPackInfoTest.cpp )
target_link_libraries ( FieldPackInfoTest walberla::blockforest walberla::core walberla::field )
waLBerla_execute_test( NAME  FieldPackInfoTest )

waLBerla_add_test_executable ( DistributionTest distributors/DistributionTest.cpp )
target_link_libraries ( DistributionTest walberla::core walberla::field )
waLBerla_execute_test( NAME DistributionTest )

waLBerla_add_test_executable ( FieldTest FieldTest.cpp )
target_link_libraries ( FieldTest walberla::core walberla::field )
waLBerla_execute_test( NAME FieldTest )

waLBerla_add_test_executable ( FieldOfCustomTypesTest FieldOfCustomTypesTest.cpp )
target_link_libraries ( FieldOfCustomTypesTest walberla::core walberla::field )
waLBerla_execute_test( NAME FieldOfCustomTypesTest )

waLBerla_add_test_executable ( FieldTiming FieldTiming.cpp )
target_link_libraries ( FieldTiming walberla::core walberla::field )
waLBerla_execute_test( NAME FieldTiming  )

waLBerla_add_test_executable ( FlagFieldTest FlagFieldTest.cpp )
target_link_libraries ( FlagFieldTest walberla::core walberla::field )
waLBerla_execute_test( NAME FlagFieldTest )

waLBerla_add_test_executable ( StabilityCheckerTest StabilityCheckerTest.cpp )
target_link_libraries ( StabilityCheckerTest walberla::blockforest walberla::field walberla::timeloop walberla::core )
waLBerla_execute_test( NAME StabilityCheckerTest )

waLBerla_add_test_executable ( InterpolationTest interpolators/InterpolationTest.cpp )
target_link_libraries ( InterpolationTest walberla::core walberla::field )
waLBerla_execute_test( NAME InterpolationTest )

waLBerla_add_test_executable ( FieldInterpolationTest interpolators/FieldInterpolationTest.cpp )
target_link_libraries ( FieldInterpolationTest walberla::core walberla::field )
waLBerla_execute_test( NAME FieldInterpolationTest )

waLBerla_add_test_executable ( AdaptorTest adaptors/AdaptorTest.cpp )
target_link_libraries ( AdaptorTest walberla::blockforest walberla::lbm walberla::core walberla::field )
waLBerla_execute_test( NAME AdaptorTest )

waLBerla_add_test_executable ( FieldGatherTest FieldGatherTest.cpp )
target_link_libraries ( FieldGatherTest walberla::blockforest walberla::core walberla::field )
waLBerla_execute_test( NAME FieldGatherTest3 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 3 DEPENDS_ON_TARGETS FieldGatherTest )
waLBerla_execute_test( NAME FieldGatherTest1 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 1 DEPENDS_ON_TARGETS FieldGatherTest )

waLBerla_add_test_executable ( FieldFileIOTest FieldFileIOTest.cpp )
target_link_libraries ( FieldFileIOTest walberla::blockforest walberla::core walberla::field )
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
    waLBerla_add_test_executable ( FieldMPIDatatypesTest FieldMPIDatatypesTest.cpp )
    target_link_libraries ( FieldMPIDatatypesTest walberla::blockforest walberla::core walberla::field )
    waLBerla_execute_test( NAME FieldMPIDatatypesTestRelease  COMMAND $<TARGET_FILE:FieldMPIDatatypesTest> PROCESSES 1                )
    waLBerla_execute_test( NAME FieldMPIDatatypesTestDebug    COMMAND $<TARGET_FILE:FieldMPIDatatypesTest> PROCESSES 1 LABELS longrun )
endif( WALBERLA_BUILD_WITH_MPI )



# CodeGen Tests
if( WALBERLA_BUILD_WITH_CODEGEN )
    waLBerla_generate_target_from_python(NAME CodegenJacobiCPUGeneratedJacobiKernel FILE codegen/JacobiKernel.py
            OUT_FILES JacobiKernel2D.cpp JacobiKernel2D.h JacobiKernel3D.cpp JacobiKernel3D.h )
    waLBerla_add_test_executable ( CodegenJacobiCPU codegen/CodegenJacobiCPU.cpp )
    target_link_libraries ( CodegenJacobiCPU walberla::timeloop CodegenJacobiCPUGeneratedJacobiKernel walberla::core
            walberla::field)
    waLBerla_execute_test( NAME CodegenJacobiCPU )

    waLBerla_generate_target_from_python(NAME SweepCollectionKernel FILE codegen/SweepCollection.py
            OUT_FILES SweepCollection.h SweepCollection.cpp)
    waLBerla_add_test_executable ( SweepCollection codegen/SweepCollection.cpp )
    target_link_libraries ( SweepCollection walberla::timeloop SweepCollectionKernel walberla::core walberla::field )
    waLBerla_execute_test( NAME SweepCollection )

    waLBerla_generate_target_from_python(NAME CodegenPoissonCPUGeneratedKernel FILE codegen/Poisson.py
            OUT_FILES Poisson.cpp Poisson.h )
    waLBerla_add_test_executable ( CodegenPoissonCPU codegen/CodegenPoissonCPU.cpp )
    target_link_libraries ( CodegenPoissonCPU walberla::timeloop CodegenPoissonCPUGeneratedKernel walberla::core
            walberla::field )
    waLBerla_execute_test( NAME CodegenPoissonCPU )

    waLBerla_generate_target_from_python(NAME CodeGenMultipleFieldSwaps FILE codegen/MultipleFieldSwaps.py
            OUT_FILES MultipleFieldSwaps.cpp MultipleFieldSwaps.h )
    waLBerla_add_test_executable ( MultipleFieldSwaps codegen/MultipleFieldSwaps.cpp )
    target_link_libraries ( MultipleFieldSwaps walberla::timeloop CodeGenMultipleFieldSwaps walberla::core
            walberla::field)
    waLBerla_execute_test( NAME MultipleFieldSwaps )

    waLBerla_generate_target_from_python(NAME CodegenGeneratedCPUFieldPackInfo FILE codegen/GeneratedFieldPackInfoTest.py
            OUT_FILES ScalarFieldCommunication.cpp ScalarFieldCommunication.h
            ScalarFieldPullReduction.cpp ScalarFieldPullReduction.h )
    waLBerla_add_test_executable ( GeneratedFieldPackInfoTest codegen/GeneratedFieldPackInfoTest.cpp )
    target_link_libraries ( GeneratedFieldPackInfoTest walberla::blockforest walberla::core walberla::field
            CodegenGeneratedCPUFieldPackInfo )
    waLBerla_execute_test( NAME GeneratedFieldPackInfoTest )

    waLBerla_generate_target_from_python(NAME CodeGenEK FILE codegen/EK.py
            OUT_FILES EKFlux.cpp EKFlux.h EKContinuity.cpp EKContinuity.h )
    waLBerla_add_test_executable ( EK codegen/EK.cpp )
    target_link_libraries ( EK walberla::blockforest walberla::core walberla::field walberla::timeloop CodeGenEK )
    waLBerla_execute_test( NAME EK )

    waLBerla_generate_target_from_python(NAME CodegenGeneratedNonUniformPackInfoTest
            FILE codegen/GeneratedFieldNonUniformPackInfoTestCPU.py
            OUT_FILES ScalarFieldNonUniformCommunicationCPU.cpp ScalarFieldNonUniformCommunicationCPU.h
    )
    waLBerla_add_test_executable ( GeneratedFieldNonUniformPackInfoTestCPU codegen/GeneratedFieldNonUniformPackInfoTestCPU.cpp )
    target_link_libraries ( GeneratedFieldNonUniformPackInfoTestCPU walberla::blockforest walberla::core walberla::field
            walberla::timeloop CodegenGeneratedNonUniformPackInfoTest )
    waLBerla_execute_test( NAME GeneratedFieldNonUniformPackInfoTestCPU )
endif()
