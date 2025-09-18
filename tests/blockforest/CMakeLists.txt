###################################################################################################
#
# Tests for walberla::blockforest  module
#
###################################################################################################

waLBerla_add_test_executable( BlockIDTest BlockIDTest.cpp )
target_link_libraries( BlockIDTest walberla::core walberla::blockforest )
waLBerla_execute_test ( NAME BlockIDTest LABELS longrun )

waLBerla_add_test_executable( SetupBlockForestTest SetupBlockForestTest.cpp )
target_link_libraries( SetupBlockForestTest walberla::core walberla::blockforest )
waLBerla_execute_test( NAME SetupBlockForestTest LABELS longrun )

waLBerla_add_test_executable( BlockForestTest BlockForestTest.cpp )
target_link_libraries( BlockForestTest walberla::core walberla::blockforest walberla::domain_decomposition )
waLBerla_execute_test( NAME BlockForestTest PROCESSES 4 )

waLBerla_add_test_executable( DeterministicCreation DeterministicCreation.cpp )
target_link_libraries( DeterministicCreation walberla::core walberla::blockforest )
waLBerla_execute_test( NAME DeterministicCreation PROCESSES 8 )

waLBerla_add_test_executable( SaveLoad SaveLoadTest.cpp )
target_link_libraries( SaveLoad walberla::core walberla::blockforest )
waLBerla_execute_test( NAME   SaveLoad01 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 1 )
waLBerla_execute_test( NAME   SaveLoad02 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 2 )
waLBerla_execute_test( NAME   SaveLoad04 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 4 )
waLBerla_execute_test( NAME   SaveLoad08 COMMAND $<TARGET_FILE:SaveLoad> PROCESSES 8 )
if( WALBERLA_BUILD_WITH_MPI )
set_property( TEST SaveLoad02 PROPERTY DEPENDS SaveLoad01 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST SaveLoad04 PROPERTY DEPENDS SaveLoad02 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
set_property( TEST SaveLoad08 PROPERTY DEPENDS SaveLoad04 ) #serialize runs of tets to avoid i/o conflicts when running ctest with -jN
endif()

waLBerla_add_test_executable( StructuredBlockForestTest StructuredBlockForestTest.cpp )
target_link_libraries( StructuredBlockForestTest walberla::core walberla::blockforest
        walberla::domain_decomposition )
waLBerla_execute_test( NAME StructuredBlockForestTest )

waLBerla_add_test_executable( BlockDataIOTest BlockDataIOTest.cpp )
target_link_libraries( BlockDataIOTest walberla::core walberla::blockforest walberla::field )
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

waLBerla_add_test_executable( GhostLayerCommTest communication/GhostLayerCommTest.cpp )
target_link_libraries( GhostLayerCommTest walberla::core walberla::blockforest walberla::field walberla::timeloop )
waLBerla_execute_test( NAME GhostLayerCommTest1 COMMAND $<TARGET_FILE:GhostLayerCommTest> DEPENDS_ON_TARGETS GhostLayerCommTest )
waLBerla_execute_test( NAME GhostLayerCommTest4 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 4 DEPENDS_ON_TARGETS GhostLayerCommTest )
waLBerla_execute_test( NAME GhostLayerCommTest8 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 8 DEPENDS_ON_TARGETS GhostLayerCommTest )

waLBerla_add_test_executable( DirectionBasedReduceCommTest communication/DirectionBasedReduceCommTest.cpp )
target_link_libraries( DirectionBasedReduceCommTest walberla::core walberla::blockforest walberla::field walberla::timeloop )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest1 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest3 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 3 DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest8 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 8 DEPENDS_ON_TARGETS DirectionBasedReduceCommTest )
