###################################################################################################
#
# Tests for walberla::domain_decomposition module
#
###################################################################################################

waLBerla_add_test_executable ( MapPointToPeriodicDomain MapPointToPeriodicDomain.cpp )
target_link_libraries( MapPointToPeriodicDomain walberla::domain_decomposition walberla::core )
waLBerla_execute_test( NAME MapPointToPeriodicDomain )

waLBerla_add_test_executable ( PeriodicIntersect PeriodicIntersect.cpp )
target_link_libraries ( PeriodicIntersect walberla::core walberla::blockforest walberla::domain_decomposition )
waLBerla_execute_test( NAME PeriodicIntersect )

waLBerla_add_test_executable ( PeriodicIntersectionVolume PeriodicIntersectionVolume.cpp )
target_link_libraries ( PeriodicIntersectionVolume walberla::domain_decomposition walberla::core walberla::blockforest )
waLBerla_execute_test( NAME PeriodicIntersectionVolume )
