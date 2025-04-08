###################################################################################################
#
# Tests for walberla::timeloop module
#
###################################################################################################

waLBerla_compile_test( FILES MultipleSweepFailTest.cpp DEPENDS walberla::blockforest )
waLBerla_execute_test( NAME MultipleSweepFailTest )
set_property( TEST MultipleSweepFailTest
      PROPERTY WILL_FAIL TRUE )

waLBerla_compile_test( FILES MultipleSweepTest.cpp DEPENDS walberla::blockforest )
waLBerla_execute_test( NAME MultipleSweepTest )

waLBerla_compile_test( FILES TimeloopSweepManagementTest.cpp DEPENDS walberla::blockforest )
waLBerla_execute_test( NAME TimeloopSweepManagementTest )

waLBerla_compile_test( FILES TimeloopAndSweepRegisterTest.cpp DEPENDS walberla::blockforest walberla::field )
waLBerla_execute_test( NAME TimeloopAndSweepRegisterTest )