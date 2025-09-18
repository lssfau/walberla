###################################################################################################
#
# Tests for walberla::core module
#
###################################################################################################

########
# cell #
########

waLBerla_add_test_executable( CellTest cell/CellTest.cpp )
target_link_libraries ( CellTest walberla::core )
waLBerla_execute_test( NAME CellTest )

waLBerla_add_test_executable( CellIntervalTest cell/CellIntervalTest.cpp )
target_link_libraries ( CellIntervalTest walberla::core )
waLBerla_execute_test( NAME CellIntervalTest )

##########
# config #
##########

waLBerla_add_test_executable( ConfigTest config/ConfigTest.cpp )
target_link_libraries ( ConfigTest walberla::core )
waLBerla_execute_test( NAME ConfigTest COMMAND $<TARGET_FILE:ConfigTest> ${CMAKE_CURRENT_SOURCE_DIR}/config/ConfigTest.dat )

#########
# debug #
#########

waLBerla_add_test_executable( DebugTest debug/DebugTest.cpp )
target_link_libraries ( DebugTest walberla::core )
waLBerla_execute_test( NAME DebugTest )

waLBerla_add_test_executable( CheckMacroTest debug/CheckMacroTest.cpp )
target_link_libraries ( CheckMacroTest walberla::core )
waLBerla_execute_test( NAME CheckMacroTest )

waLBerla_add_test_executable( PrintStacktraceTest debug/PrintStacktraceTest.cpp )
target_link_libraries ( PrintStacktraceTest walberla::core )

###########
# logging #
###########

waLBerla_add_test_executable( LoggingTest logging/LoggingTest.cpp )
target_link_libraries ( LoggingTest walberla::core )
waLBerla_execute_test( NAME LoggingTest )

########
# math #
########

waLBerla_add_test_executable( PrimesTest math/PrimesTest.cpp )
target_link_libraries ( PrimesTest walberla::core )
waLBerla_execute_test( NAME PrimesTest )

waLBerla_add_test_executable( RandomTest math/RandomTest.cpp )
target_link_libraries ( RandomTest walberla::core )
waLBerla_execute_test( NAME RandomTest )

waLBerla_add_test_executable( SampleTest math/SampleTest.cpp )
target_link_libraries ( SampleTest walberla::core )
waLBerla_execute_test( NAME SampleTest1Proc  COMMAND $<TARGET_FILE:SampleTest> )
waLBerla_execute_test( NAME SampleTest4Proc  COMMAND $<TARGET_FILE:SampleTest> PROCESSES 4  )
waLBerla_execute_test( NAME SampleTest10Proc COMMAND $<TARGET_FILE:SampleTest> PROCESSES 10 LABELS longrun )

waLBerla_add_test_executable( DistributedSampleTest math/DistributedSampleTest.cpp )
target_link_libraries ( DistributedSampleTest walberla::core )
waLBerla_execute_test( NAME DistributedSampleTest1Proc  COMMAND $<TARGET_FILE:DistributedSampleTest> )
waLBerla_execute_test( NAME DistributedSampleTest4Proc  COMMAND $<TARGET_FILE:DistributedSampleTest> PROCESSES 4  )
waLBerla_execute_test( NAME DistributedSampleTest10Proc COMMAND $<TARGET_FILE:DistributedSampleTest> PROCESSES 10 LABELS longrun )

waLBerla_add_test_executable( ParserTest math/ParserTest.cpp )
target_link_libraries ( ParserTest walberla::core )
waLBerla_execute_test( NAME ParserTest )

waLBerla_add_test_executable( PlaneTest math/PlaneTest.cpp )
target_link_libraries ( PlaneTest walberla::core )
waLBerla_execute_test( NAME PlaneTest )

waLBerla_add_test_executable( MatrixVector2Test math/MatrixVector2Test.cpp )
target_link_libraries ( MatrixVector2Test walberla::core )
waLBerla_execute_test( NAME MatrixVector2Test )


waLBerla_add_test_executable( Matrix3Test math/Matrix3Test.cpp )
target_link_libraries ( Matrix3Test walberla::core )
waLBerla_execute_test( NAME Matrix3Test )

waLBerla_add_test_executable( GenericAABBTest math/GenericAABBTest.cpp )
target_link_libraries ( GenericAABBTest walberla::core walberla::stencil walberla::domain_decomposition )
waLBerla_execute_test( NAME GenericAABBTest )


waLBerla_add_test_executable( FastInvSqrtTest math/FastInvSqrtTest.cpp )
target_link_libraries ( FastInvSqrtTest walberla::core )
waLBerla_execute_test( NAME FastInvSqrtTest )

#######
# mpi #
#######

waLBerla_add_test_executable( BufferTest mpi/BufferTest.cpp )
target_link_libraries ( BufferTest walberla::core )
waLBerla_execute_test( NAME BufferTest )

waLBerla_add_test_executable( BufferIOTest mpi/BufferIOTest.cpp )
target_link_libraries ( BufferIOTest walberla::core )
waLBerla_execute_test( NAME BufferIOTest1 COMMAND $<TARGET_FILE:BufferIOTest> )
waLBerla_execute_test( NAME BufferIOTest4 COMMAND $<TARGET_FILE:BufferIOTest> PROCESSES 4)
waLBerla_execute_test( NAME BufferIOTest8 COMMAND $<TARGET_FILE:BufferIOTest> PROCESSES 8)
if( WALBERLA_BUILD_WITH_MPI )
set_property( TEST BufferIOTest4 PROPERTY DEPENDS BufferIOTest1 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST BufferIOTest8 PROPERTY DEPENDS BufferIOTest4 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
endif()

waLBerla_add_test_executable( BufferSystemTest mpi/BufferSystemTest.cpp )
target_link_libraries ( BufferSystemTest walberla::core )
waLBerla_execute_test( NAME BufferSystemTest PROCESSES 4 )

waLBerla_add_test_executable( BroadcastTest mpi/BroadcastTest.cpp )
target_link_libraries ( BroadcastTest walberla::core )
waLBerla_execute_test( NAME BroadcastTest1 COMMAND $<TARGET_FILE:BroadcastTest> )
waLBerla_execute_test( NAME BroadcastTest4 COMMAND $<TARGET_FILE:BroadcastTest> PROCESSES 4)

waLBerla_add_test_executable( ReduceTest mpi/ReduceTest.cpp )
target_link_libraries ( ReduceTest walberla::core )
waLBerla_execute_test( NAME ReduceTest1 COMMAND $<TARGET_FILE:ReduceTest> )
waLBerla_execute_test( NAME ReduceTest4 COMMAND $<TARGET_FILE:ReduceTest> PROCESSES 4)

waLBerla_add_test_executable( GatherTest mpi/GatherTest.cpp )
target_link_libraries ( GatherTest walberla::core walberla::gather )
waLBerla_execute_test( NAME GatherTest1 COMMAND $<TARGET_FILE:GatherTest> )
waLBerla_execute_test( NAME GatherTest4 COMMAND $<TARGET_FILE:GatherTest> PROCESSES 4)

waLBerla_add_test_executable( GathervTest mpi/GathervTest.cpp )
target_link_libraries ( GathervTest walberla::core walberla::gather )
waLBerla_execute_test( NAME GathervTest1 COMMAND $<TARGET_FILE:GathervTest> )
waLBerla_execute_test( NAME GathervTest4 COMMAND $<TARGET_FILE:GathervTest> PROCESSES 4)

waLBerla_add_test_executable( MPITextFileTest mpi/MPITextFileTest.cpp )
target_link_libraries ( MPITextFileTest walberla::core )
waLBerla_execute_test( NAME MPITextFileTest1 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_1.txt 16 )
waLBerla_execute_test( NAME MPITextFileTest4 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_4.txt 16 PROCESSES 4 )
waLBerla_execute_test( NAME MPITextFileTest8 COMMAND $<TARGET_FILE:MPITextFileTest> MPI_Testfile_8.txt 16 PROCESSES 8 )

waLBerla_add_test_executable( SetReductionTest mpi/SetReductionTest.cpp )
target_link_libraries ( SetReductionTest walberla::core )
waLBerla_execute_test( NAME SetReductionTest1  COMMAND $<TARGET_FILE:SetReductionTest> )
waLBerla_execute_test( NAME SetReductionTest2  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 2 )
waLBerla_execute_test( NAME SetReductionTest3  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 3 )
waLBerla_execute_test( NAME SetReductionTest4  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 4 )
waLBerla_execute_test( NAME SetReductionTest5  COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 5 )
waLBerla_execute_test( NAME SetReductionTest27 COMMAND $<TARGET_FILE:SetReductionTest> PROCESSES 27 )

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
   waLBerla_add_test_executable( MPIFloat16Test mpi/MPIFloat16Test.cpp )
   target_link_libraries ( MPIFloat16Test walberla::core )
   target_compile_features( MPIFloat16Test PUBLIC cxx_std_23 )
   waLBerla_execute_test( NAME MPIFloat16Test4 COMMAND $<TARGET_FILE:MPIFloat16Test> PROCESSES 4)
endif ()


##############
# selectable #
##############

waLBerla_add_test_executable( SetSelectableObjectTest selectable/SetSelectableObjectTest.cpp )
target_link_libraries ( SetSelectableObjectTest walberla::core )
waLBerla_execute_test( NAME SetSelectableObjectTest  )

##########
# timing #
##########

waLBerla_add_test_executable( ParallelTimerTest timing/ParallelTimerTest.cpp )
target_link_libraries ( ParallelTimerTest walberla::core )
waLBerla_execute_test( NAME ParallelTimerTest PROCESSES 2 )

waLBerla_add_test_executable( TimerTest timing/TimerTest.cpp )
target_link_libraries ( TimerTest walberla::core )
waLBerla_execute_test( NAME TimerTest )

waLBerla_add_test_executable( TimingPoolTest timing/TimingPoolTest.cpp )
target_link_libraries ( TimingPoolTest walberla::core )
waLBerla_execute_test( NAME TimingPoolTest PROCESSES 3  )

waLBerla_add_test_executable( TimingTreeTest timing/TimingTreeTest.cpp )
target_link_libraries ( TimingTreeTest walberla::core )
waLBerla_execute_test( NAME TimingTreeTest PROCESSES 3)

waLBerla_add_test_executable( SweepTimeloopTimerReduction timing/SweepTimeloopTimerReduction.cpp )
target_link_libraries ( SweepTimeloopTimerReduction walberla::core walberla::blockforest walberla::domain_decomposition
        walberla::timeloop )
waLBerla_execute_test( NAME SweepTimeloopTimerReduction PROCESSES 9  )


########
# walberla::core #
########

if ( WALBERLA_BUILD_WITH_PYTHON )
   waLBerla_add_test_executable( AllHeaderTest AllHeaderTest.cpp )
   target_link_libraries ( AllHeaderTest walberla::core walberla::blockforest walberla::field walberla::geometry
           walberla::python_coupling )
else()
   waLBerla_add_test_executable( AllHeaderTest AllHeaderTest.cpp )
   target_link_libraries ( AllHeaderTest walberla::core walberla::blockforest walberla::field walberla::geometry )
endif()
waLBerla_execute_test( NAME AllHeaderTest )

waLBerla_add_test_executable( Angles Angles.test.cpp )
target_link_libraries ( Angles walberla::core )
waLBerla_execute_test( NAME Angles )

waLBerla_add_test_executable( ConcatIterator ConcatIterator.cpp )
target_link_libraries ( ConcatIterator walberla::core )
waLBerla_execute_test( NAME ConcatIterator )

waLBerla_add_test_executable( DataTypesTest DataTypesTest.cpp )
target_link_libraries ( DataTypesTest walberla::core )
waLBerla_execute_test( NAME DataTypesTest )

waLBerla_add_test_executable( DebugSTLTest DebugSTLTest.cpp )
target_link_libraries ( DebugSTLTest walberla::core )
waLBerla_execute_test( NAME DebugSTLTest )
set_tests_properties(DebugSTLTest PROPERTIES WILL_FAIL TRUE)

waLBerla_add_test_executable( FunctionTraitsTest FunctionTraitsTest.cpp )
target_link_libraries ( FunctionTraitsTest walberla::core )
waLBerla_execute_test( NAME FunctionTraitsTest )

waLBerla_add_test_executable( GridGeneratorTest GridGeneratorTest.cpp )
target_link_libraries ( GridGeneratorTest walberla::core )
waLBerla_execute_test( NAME GridGeneratorTest )

waLBerla_add_test_executable( MemoryUsage MemoryUsage.cpp )
target_link_libraries ( MemoryUsage walberla::core )
waLBerla_execute_test( NAME MemoryUsage )

waLBerla_add_test_executable( OpenMPWrapperTest OpenMPWrapperTest.cpp )
target_link_libraries ( OpenMPWrapperTest walberla::core )
waLBerla_execute_test( NAME OpenMPWrapperTest )

waLBerla_add_test_executable( RandomUUID RandomUUID.cpp )
target_link_libraries ( RandomUUID walberla::core )
waLBerla_execute_test( NAME RandomUUID )

waLBerla_add_test_executable( SetTest SetTest.cpp )
target_link_libraries ( SetTest walberla::core )
waLBerla_execute_test( NAME SetTest )

waLBerla_add_test_executable( UNIQUEID UniqueID.cpp )
target_link_libraries ( UNIQUEID walberla::core )
waLBerla_execute_test( NAME UNIQUEID PROCESSES 4)

waLBerla_add_test_executable( VersionTest VersionTest.cpp )
target_link_libraries ( VersionTest walberla::core )
waLBerla_execute_test( NAME VersionTest )

if ( WALBERLA_CXX_COMPILER_IS_INTEL AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "20.0" )
   waLBerla_add_test_executable( VariantTest VariantTest.cpp )
   target_link_libraries ( VariantTest walberla::core )
   waLBerla_execute_test( NAME VariantTest )
endif()

##################
# load_balancing #
##################

if( WALBERLA_BUILD_WITH_METIS )
   waLBerla_add_test_executable( MetisTest load_balancing/MetisTest.cpp )
   target_link_libraries ( MetisTest walberla::core walberla::field )
   waLBerla_execute_test( NAME MetisTest COMMAND $<TARGET_FILE:MetisTest> 64 64  4 --no-vtk )
endif()

if( WALBERLA_BUILD_WITH_PARMETIS )
   waLBerla_add_test_executable( ParMetisTest load_balancing/ParMetisTest.cpp )
   target_link_libraries ( ParMetisTest walberla::core walberla::blockforest walberla::field walberla::stencil
           walberla::vtk )
   waLBerla_execute_test( NAME ParMetisTest1 COMMAND $<TARGET_FILE:ParMetisTest> 64 64  4 --no-vtk )
   waLBerla_execute_test( NAME ParMetisTest2 COMMAND $<TARGET_FILE:ParMetisTest> 64 64  8 --no-vtk PROCESSES 2 )
   waLBerla_execute_test( NAME ParMetisTest4 COMMAND $<TARGET_FILE:ParMetisTest> 64 64 16 --no-vtk PROCESSES 4 )

   waLBerla_add_test_executable( PlainParMetisTest load_balancing/PlainParMetisTest.cpp )
   target_link_libraries ( PlainParMetisTest walberla::core )
   waLBerla_execute_test( NAME PlainParMetisTest PROCESSES 3 )
endif()

###################
# Mixed Precision #
###################

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
   waLBerla_add_test_executable( Float16SupportTest Float16SupportTest.cpp )
   target_link_libraries ( Float16SupportTest walberla::core )
   # Actual support for float16 is available only since C++23
   #   before is_arithmetic and is_floating_point evaluated to false,
   #   also many STL functions are compatible with float16 only since C++23.
   # Which features are actually supported depend on the compiler
   target_compile_features( Float16SupportTest PUBLIC cxx_std_23 )
   waLBerla_execute_test(NAME Float16SupportTest)

   waLBerla_add_test_executable( FP16Test FP16Test.cpp )
   target_link_libraries ( FP16Test walberla::core )
   waLBerla_execute_test( NAME FP16Test )
endif ()