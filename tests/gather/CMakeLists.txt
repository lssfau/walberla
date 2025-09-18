###################################################################################################
#
# Tests for gather module
#
###################################################################################################


waLBerla_add_test_executable ( MPIGatherSchemeTest MPIGatherSchemeTest.cpp )
target_link_libraries ( MPIGatherSchemeTest walberla::core walberla::gather walberla::blockforest )
waLBerla_execute_test( NAME MPIGatherSchemeTest PROCESSES 7 )

waLBerla_add_test_executable ( CurveGatherTest CurveGatherTest.cpp )
target_link_libraries ( CurveGatherTest walberla::core walberla::gather walberla::blockforest walberla::field
        walberla::timeloop )
waLBerla_execute_test( NAME CurveGatherTest PROCESSES 4 )

waLBerla_add_test_executable ( GatherSchemeTest GatherSchemeTest.cpp )
target_link_libraries ( GatherSchemeTest walberla::core walberla::gather walberla::blockforest walberla::field )
waLBerla_execute_test( NAME GatherSchemeTest PROCESSES 4 )
