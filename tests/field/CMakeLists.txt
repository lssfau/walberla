###################################################################################################
#
# Tests for field module
#
###################################################################################################

waLBerla_compile_test( FILES AccuracyEvaluationTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME AccuracyEvaluationTest4 COMMAND $<TARGET_FILE:AccuracyEvaluationTest> PROCESSES 4  )

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

waLBerla_compile_test( FILES interpolators/InterpolationTest.cpp)
waLBerla_execute_test( NAME InterpolationTest )

waLBerla_compile_test( FILES interpolators/FieldInterpolationTest.cpp)
waLBerla_execute_test( NAME FieldInterpolationTest )

waLBerla_compile_test( FILES adaptors/AdaptorTest.cpp DEPENDS blockforest gui lbm )
waLBerla_execute_test( NAME AdaptorTest )

waLBerla_compile_test( FILES FieldGatherTest.cpp DEPENDS blockforest gui )
waLBerla_execute_test( NAME FieldGatherTest3 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 3 )
waLBerla_execute_test( NAME FieldGatherTest1 COMMAND $<TARGET_FILE:FieldGatherTest>  PROCESSES 1 )

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

waLBerla_python_file_generates(codegen/JacobiKernel.py JacobiKernel2D.cpp JacobiKernel3D.cpp )
waLBerla_compile_test( FILES codegen/CodegenJacobiCPU.cpp codegen/JacobiKernel.py DEPENDS gui timeloop )
waLBerla_execute_test( NAME CodegenJacobiCPU )


