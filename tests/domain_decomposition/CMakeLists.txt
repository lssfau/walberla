###################################################################################################
#
# Tests for walberla::domain_decomposition module
#
###################################################################################################

waLBerla_compile_test( NAME MapPointToPeriodicDomain FILES MapPointToPeriodicDomain.cpp DEPENDS walberla::core )
waLBerla_execute_test( NAME MapPointToPeriodicDomain )

waLBerla_compile_test( NAME PeriodicIntersect FILES PeriodicIntersect.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME PeriodicIntersect )

waLBerla_compile_test( NAME PeriodicIntersectionVolume FILES PeriodicIntersectionVolume.cpp DEPENDS walberla::core walberla::blockforest )
waLBerla_execute_test( NAME PeriodicIntersectionVolume )
