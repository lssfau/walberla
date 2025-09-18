###################################################################################################
#
# Tests for walberla::sqlite module
#
###################################################################################################

waLBerla_add_test_executable ( SQLiteTest SQLiteTest.cpp )
target_link_libraries ( SQLiteTest walberla::sqlite )
waLBerla_execute_test( NAME SQLiteTest )