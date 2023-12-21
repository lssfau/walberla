###################################################################################################
#
# Tests for core module
#
###################################################################################################

########
# cell #
########

waLBerla_compile_test( FILES cell/CellTest.cpp )
waLBerla_execute_test( NAME CellTest )

waLBerla_compile_test( FILES cell/CellIntervalTest.cpp )
waLBerla_execute_test( NAME CellIntervalTest )

##########
# config #
##########

waLBerla_compile_test( FILES config/ConfigTest.cpp )
waLBerla_execute_test( NAME ConfigTest COMMAND $<TARGET_FILE:ConfigTest> ${CMAKE_CURRENT_SOURCE_DIR}/config/ConfigTest.dat )

#########
# debug #
#########

waLBerla_compile_test( FILES debug/DebugTest.cpp)
waLBerla_execute_test( NAME DebugTest )

waLBerla_compile_test( FILES debug/CheckMacroTest.cpp )
waLBerla_execute_test( NAME CheckMacroTest )

waLBerla_compile_test( FILES debug/PrintStacktraceTest.cpp )

###########
# logging #
###########

waLBerla_compile_test( FILES logging/LoggingTest.cpp)
waLBerla_execute_test( NAME LoggingTest )

########
# math #
########

waLBerla_compile_test( FILES math/PrimesTest.cpp )
waLBerla_execute_test( NAME PrimesTest )

waLBerla_compile_test( FILES math/RandomTest.cpp )
waLBerla_execute_test( NAME RandomTest )

waLBerla_compile_test( FILES math/SampleTest.cpp )
waLBerla_execute_test( NAME SampleTest1Proc  COMMAND $<TARGET_FILE:SampleTest> )
waLBerla_execute_test( NAME SampleTest4Proc  COMMAND $<TARGET_FILE:SampleTest> PROCESSES 4  )
waLBerla_execute_test( NAME SampleTest10Proc COMMAND $<TARGET_FILE:SampleTest> PROCESSES 10 LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES math/DistributedSampleTest.cpp )
waLBerla_execute_test( NAME DistributedSampleTest1Proc  COMMAND $<TARGET_FILE:DistributedSampleTest> )
waLBerla_execute_test( NAME DistributedSampleTest4Proc  COMMAND $<TARGET_FILE:DistributedSampleTest> PROCESSES 4  )
waLBerla_execute_test( NAME DistributedSampleTest10Proc COMMAND $<TARGET_FILE:DistributedSampleTest> PROCESSES 10 LABELS longrun )

waLBerla_compile_test( FILES math/ParserTest.cpp)
waLBerla_execute_test( NAME ParserTest )

waLBerla_compile_test( FILES math/PlaneTest.cpp )
waLBerla_execute_test( NAME PlaneTest )

waLBerla_compile_test( FILES math/MatrixVector2Test.cpp )
waLBerla_execute_test( NAME MatrixVector2Test )


waLBerla_compile_test( FILES math/Matrix3Test.cpp )
waLBerla_execute_test( NAME Matrix3Test )

waLBerla_compile_test( FILES math/GenericAABBTest.cpp DEPENDS stencil domain_decomposition )
waLBerla_execute_test( NAME GenericAABBTest )


waLBerla_compile_test( FILES math/FastInvSqrtTest.cpp )
waLBerla_execute_test( NAME FastInvSqrtTest )

#######
# mpi #
#######

waLBerla_compile_test( FILES mpi/BufferTest.cpp )
waLBerla_execute_test( NAME BufferTest )

waLBerla_compile_test( FILES mpi/BufferSystemTest.cpp )
waLBerla_execute_test( NAME BufferSystemTest PROCESSES 4 )

waLBerla_compile_test( FILES mpi/BroadcastTest.cpp )
waLBerla_execute_test( NAME BroadcastTest1 COMMAND $<TARGET_FILE:BroadcastTest> )
waLBerla_execute_test( NAME BroadcastTest4 COMMAND $<TARGET_FILE:BroadcastTest> PROCESSES 4)

waLBerla_compile_test( FILES mpi/ReduceTest.cpp )
waLBerla_execute_test( NAME ReduceTest1 COMMAND $<TARGET_FILE:ReduceTest> )
waLBerla_execute_test( NAME ReduceTest4 COMMAND $<TARGET_FILE:ReduceTest> PROCESSES 4)

waLBerla_compile_test( FILES mpi/GatherTest.cpp )
waLBerla_execute_test( NAME GatherTest1 COMMAND $<TARGET_FILE:GatherTest> )
waLBerla_execute_test( NAME GatherTest4 COMMAND $<TARGET_FILE:GatherTest> PROCESSES 4)

waLBerla_compile_test( FILES mpi/GathervTest.cpp )
waLBerla_execute_test( NAME GathervTest1 COMMAND $<TARGET_FILE:GathervTest> )
waLBerla_execute_test( NAME GathervTest4 COMMAND $<TARGET_FILE:GathervTest> PROCESSES 4)

waLBerla_compile_test( FILES mpi/MPITextFileTest.cpp )
waLBerla_execute_test( NAME MPITextFileTest1 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_1.txt 16 )
waLBerla_execute_test( NAME MPITextFileTest4 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_4.txt 16 PROCESSES 4 )
waLBerla_execute_test( NAME MPITextFileTest8 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_8.txt 16 PROCESSES 8 )

waLBerla_compile_test( FILES mpi/SetReductionTest.cpp )
waLBerla_execute_test( NAME SetReductionTest1  COMMAND $<TARGET_FILE:SetReductionTest> )
waLBerla_execute_test( NAME SetReductionTest2  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 2 )
waLBerla_execute_test( NAME SetReductionTest3  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 3 )
waLBerla_execute_test( NAME SetReductionTest4  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 4 )
waLBerla_execute_test( NAME SetReductionTest5  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 5 )
waLBerla_execute_test( NAME SetReductionTest27 COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 27 )



##############
# selectable #
##############

waLBerla_compile_test( FILES selectable/SetSelectableObjectTest.cpp )
waLBerla_execute_test( NAME SetSelectableObjectTest  )

##########
# timing #
##########

waLBerla_compile_test( FILES timing/ParallelTimerTest.cpp )
waLBerla_execute_test( NAME ParallelTimerTest PROCESSES 2 )

waLBerla_compile_test( FILES timing/TimerTest.cpp )
waLBerla_execute_test( NAME TimerTest )

waLBerla_compile_test( FILES timing/TimingPoolTest.cpp )
waLBerla_execute_test( NAME TimingPoolTest PROCESSES 3  )

waLBerla_compile_test( FILES timing/TimingTreeTest.cpp )
waLBerla_execute_test( NAME TimingTreeTest PROCESSES 3)

waLBerla_compile_test( FILES timing/SweepTimeloopTimerReduction.cpp DEPENDS blockforest domain_decomposition timeloop)
waLBerla_execute_test( NAME SweepTimeloopTimerReduction PROCESSES 9  )


########
# core #
########

if ( WALBERLA_BUILD_WITH_PYTHON )
   waLBerla_compile_test( FILES AllHeaderTest.cpp DEPENDS blockforest field geometry pe pe_coupling python_coupling )
else()
   waLBerla_compile_test( FILES AllHeaderTest.cpp DEPENDS blockforest field geometry pe pe_coupling )
endif()
waLBerla_execute_test( NAME AllHeaderTest )

waLBerla_compile_test( FILES Angles.test.cpp )
waLBerla_execute_test( NAME Angles )

waLBerla_compile_test( FILES ConcatIterator.cpp )
waLBerla_execute_test( NAME ConcatIterator )

waLBerla_compile_test( FILES DataTypesTest.cpp )
waLBerla_execute_test( NAME DataTypesTest )

waLBerla_compile_test( FILES DebugSTLTest.cpp )
waLBerla_execute_test( NAME DebugSTLTest )
set_tests_properties(DebugSTLTest PROPERTIES WILL_FAIL TRUE)

waLBerla_compile_test( FILES FP16Test.cpp )
waLBerla_execute_test( NAME FP16Test )

waLBerla_compile_test( FILES FunctionTraitsTest.cpp )
waLBerla_execute_test( NAME FunctionTraitsTest )

waLBerla_compile_test( FILES GridGeneratorTest.cpp )
waLBerla_execute_test( NAME GridGeneratorTest )

waLBerla_compile_test( FILES MemoryUsage.cpp )
waLBerla_execute_test( NAME MemoryUsage )

waLBerla_compile_test( FILES OpenMPWrapperTest.cpp )
waLBerla_execute_test( NAME OpenMPWrapperTest )

waLBerla_compile_test( FILES RandomUUID.cpp )
waLBerla_execute_test( NAME RandomUUID )

waLBerla_compile_test( FILES SetTest.cpp )
waLBerla_execute_test( NAME SetTest )

waLBerla_compile_test( NAME UNIQUEID FILES UniqueID.cpp )
waLBerla_execute_test( NAME UNIQUEID PROCESSES 4)

waLBerla_compile_test( FILES VersionTest.cpp )
waLBerla_execute_test( NAME VersionTest )

if ( WALBERLA_CXX_COMPILER_IS_INTEL AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "20.0" )
  waLBerla_compile_test( FILES VariantTest.cpp )
  waLBerla_execute_test( NAME VariantTest )
endif()

##################
# load_balancing #
##################

if( WALBERLA_BUILD_WITH_METIS )
   waLBerla_compile_test( NAME MetisTest FILES load_balancing/MetisTest.cpp DEPENDS field )
   waLBerla_execute_test( NAME MetisTest COMMAND $<TARGET_FILE:MetisTest> 64 64  4 --no-vtk )
endif()

if( WALBERLA_BUILD_WITH_PARMETIS )
   waLBerla_compile_test( NAME ParMetisTest FILES load_balancing/ParMetisTest.cpp DEPENDS blockforest field stencil vtk  )
   waLBerla_execute_test( NAME ParMetisTest1 COMMAND $<TARGET_FILE:ParMetisTest> 64 64  4 --no-vtk )
   waLBerla_execute_test( NAME ParMetisTest2 COMMAND $<TARGET_FILE:ParMetisTest> 64 64  8 --no-vtk PROCESSES 2 )
   waLBerla_execute_test( NAME ParMetisTest4 COMMAND $<TARGET_FILE:ParMetisTest> 64 64 16 --no-vtk PROCESSES 4 )

   waLBerla_compile_test( NAME PlainParMetisTest FILES load_balancing/PlainParMetisTest.cpp )
   waLBerla_execute_test( NAME PlainParMetisTest PROCESSES 3 )
endif()

###################
# Mixed Precision #
###################

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
   waLBerla_compile_test( Name Float16SupportTest FILES Float16SupportTest.cpp DEPENDS core)
   # Actual support for float16 is available only since C++23
   #   before is_arithmetic and is_floating_point evaluated to false,
   #   also many STL functions are compatible with float16 only since C++23.
   # Which features are actually supported depend on the compiler
   target_compile_features( Float16SupportTest PUBLIC cxx_std_23 )
   waLBerla_execute_test(NAME Float16SupportTest)
endif ()