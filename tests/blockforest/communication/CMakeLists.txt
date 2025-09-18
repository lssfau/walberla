###################################################################################################
#
# Tests for comm module
#
###################################################################################################

waLBerla_add_test_executable ( GhostLayerCommTest GhostLayerCommTest.cpp )
target_link_libraries ( GhostLayerCommTest walberla::timeloop walberla::communication )
waLBerla_execute_test( NAME GhostLayerCommTest1 COMMAND $<TARGET_FILE:GhostLayerCommTest> )
waLBerla_execute_test( NAME GhostLayerCommTest4 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 4 )
waLBerla_execute_test( NAME GhostLayerCommTest8 COMMAND $<TARGET_FILE:GhostLayerCommTest> PROCESSES 8 )

waLBerla_add_test_executable ( DirectionBasedReduceCommTest DirectionBasedReduceCommTest.cpp )
target_link_libraries ( DirectionBasedReduceCommTest walberla::timeloop walberla::communication )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest1 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest3 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 3 )
waLBerla_execute_test( NAME DirectionBasedReduceCommTest8 COMMAND $<TARGET_FILE:DirectionBasedReduceCommTest> PROCESSES 8 )
