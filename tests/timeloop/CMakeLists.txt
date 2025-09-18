###################################################################################################
#
# Tests for walberla::timeloop module
#
###################################################################################################

waLBerla_add_test_executable ( MultipleSweepFailTest MultipleSweepFailTest.cpp )
target_link_libraries ( MultipleSweepFailTest walberla::blockforest walberla::timeloop )
waLBerla_execute_test( NAME MultipleSweepFailTest )
set_property( TEST MultipleSweepFailTest
      PROPERTY WILL_FAIL TRUE )

waLBerla_add_test_executable ( MultipleSweepTest MultipleSweepTest.cpp )
target_link_libraries ( MultipleSweepTest walberla::blockforest walberla::timeloop )
waLBerla_execute_test( NAME MultipleSweepTest )

waLBerla_add_test_executable ( TimeloopSweepManagementTest TimeloopSweepManagementTest.cpp )
target_link_libraries ( TimeloopSweepManagementTest walberla::blockforest walberla::timeloop )
waLBerla_execute_test( NAME TimeloopSweepManagementTest )

waLBerla_add_test_executable ( TimeloopAndSweepRegisterTest TimeloopAndSweepRegisterTest.cpp )
target_link_libraries ( TimeloopAndSweepRegisterTest walberla::blockforest walberla::timeloop walberla::field )
waLBerla_execute_test( NAME TimeloopAndSweepRegisterTest )