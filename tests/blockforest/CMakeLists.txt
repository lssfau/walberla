###################################################################################################
#
# Tests for blockforest module
#
###################################################################################################

waLBerla_compile_test ( FILES BlockIDTest.cpp )
waLBerla_execute_test ( NAME BlockIDTest LABELS longrun )

waLBerla_compile_test( FILES SetupBlockForestTest.cpp )
waLBerla_execute_test( NAME SetupBlockForestTest LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES BlockForestTest.cpp )
waLBerla_execute_test( NAME BlockForestTest PROCESSES 4 )

waLBerla_compile_test( NAME   SaveLoad FILES SaveLoadTest.cpp DEPENDS core blockforest  )
waLBerla_execute_test( NAME   SaveLoad01 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 1 )
waLBerla_execute_test( NAME   SaveLoad02 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 2 )
waLBerla_execute_test( NAME   SaveLoad04 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 4 )
waLBerla_execute_test( NAME   SaveLoad08 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 8 )
if( WALBERLA_BUILD_WITH_MPI )
set_property( TEST SaveLoad02 PROPERTY DEPENDS SaveLoad01 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST SaveLoad04 PROPERTY DEPENDS SaveLoad02 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST SaveLoad08 PROPERTY DEPENDS SaveLoad04 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
endif()

waLBerla_compile_test( FILES StructuredBlockForestTest.cpp )
waLBerla_execute_test( NAME StructuredBlockForestTest )

waLBerla_compile_test( FILES BlockDataIOTest.cpp DEPENDS field )
waLBerla_execute_test( NAME BlockDataIOTest1 COMMAND $<TARGET_FILE:BlockDataIOTest> )
waLBerla_execute_test( NAME BlockDataIOTest3 COMMAND $<TARGET_FILE:BlockDataIOTest> PROCESSES 3 )
waLBerla_execute_test( NAME BlockDataIOTest8 COMMAND $<TARGET_FILE:BlockDataIOTest> PROCESSES 8 )
#serialize runs of tests to avoid i/o conflicts when running ctest with -jN
if( WALBERLA_BUILD_WITH_MPI )
   set_property( TEST BlockDataIOTest3 PROPERTY DEPENDS BlockDataIOTest1 )
   set_property( TEST BlockDataIOTest8 PROPERTY DEPENDS BlockDataIOTest3 )
endif( WALBERLA_BUILD_WITH_MPI )

# communication

waLBerla_compile_test( FILES communication/GhostLayerCommTest.cpp DEPENDS field timeloop )
waLBerla_execute_test( NAME GhostLayerCommTest1 COMMAND $<TARGET_FILE:GhostLayerCommTest> )
waLBerla_execute_test( NAME GhostLayerCommTest4 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 4 )
waLBerla_execute_test( NAME GhostLayerCommTest8 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 8 )

waLBerla_compile_test( FILES communication/DirectionBasedReduceCommTest.cpp DEPENDS field timeloop )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest1 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest3 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 3 )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest8 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 8 )
