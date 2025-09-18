###################################################################################################
#
# Tests for executiontree module
#
###################################################################################################

waLBerla_add_test_executable( ExecutionTreeTest ExecutionTreeTest.cpp )
target_link_libraries ( ExecutionTreeTest walberla::executiontree )
waLBerla_execute_test( NAME ExecutionTreeTest )
