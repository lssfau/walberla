###################################################################################################
#
# Tests for blockforest module
#
###################################################################################################

waLBerla_compile_test ( FILES BlockIDTest.cpp DEPENDS core )
waLBerla_execute_test ( NAME BlockIDTest LABELS longrun )

waLBerla_compile_test( FILES SetupBlockForestTest.cpp DEPENDS core )
waLBerla_execute_test( NAME SetupBlockForestTest LABELS longrun CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_compile_test( FILES BlockForestTest.cpp DEPENDS core domain_decomposition )
waLBerla_execute_test( NAME BlockForestTest PROCESSES 4 )

waLBerla_compile_test( FILES DeterministicCreation.cpp DEPENDS core )
waLBerla_execute_test( NAME DeterministicCreation PROCESSES 8 )

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

waLBerla_compile_test( FILES StructuredBlockForestTest.cpp DEPENDS core domain_decomposition )
waLBerla_execute_test( NAME StructuredBlockForestTest )

waLBerla_compile_test( FILES BlockDataIOTest.cpp DEPENDS core field )
waLBerla_execute_test( NAME BlockDataIOTest1 COMMAND $<TARGET_FILE:BlockDataIOTest> DEPENDS_ON_TARGETS BlockDataIOTest )
waLBerla_execute_test( NAME BlockDataIOTest3 COMMAND $<TARGET_FILE:BlockDataIOTest> PROCESSES 3 DEPENDS_ON_TARGETS BlockDataIOTest )
waLBerla_execute_test( NAME BlockDataIOTest8 COMMAND $<TARGET_FILE:BlockDataIOTest> PROCESSES 8 DEPENDS_ON_TARGETS BlockDataIOTest )
#serialize runs of tests to avoid i/o conflicts when running ctest with -jN
if( WALBERLA_BUILD_WITH_MPI )
   if( TARGET BlockDataIOTest )
      set_property( TEST BlockDataIOTest3 PROPERTY DEPENDS BlockDataIOTest1 )
      set_property( TEST BlockDataIOTest8 PROPERTY DEPENDS BlockDataIOTest3 )
   endif()
endif( WALBERLA_BUILD_WITH_MPI )

# communication

waLBerla_compile_test( FILES communication/GhostLayerCommTest.cpp DEPENDS field timeloop )
waLBerla_execute_test( NAME GhostLayerCommTest1 COMMAND $<TARGET_FILE:GhostLayerCommTest> DEPENDS_ON_TARGETS GhostLayerCommTest )
waLBerla_execute_test( NAME GhostLayerCommTest4 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 4 DEPENDS_ON_TARGETS GhostLayerCommTest )
waLBerla_execute_test( NAME GhostLayerCommTest8 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 8 DEPENDS_ON_TARGETS GhostLayerCommTest )

waLBerla_compile_test( FILES communication/DirectionBasedReduceCommTest.cpp DEPENDS field timeloop )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest1 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest3 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 3 DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest8 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 8 DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
