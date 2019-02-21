###################################################################################################
#
# Tests for executiontree module
#
###################################################################################################

waLBerla_compile_test( NAME ExecutionTreeTest FILES ExecutionTreeTest.cpp )
waLBerla_execute_test( NAME ExecutionTreeTest )
