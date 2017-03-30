###################################################################################################
#
# Tests for domain_decomposition module
#
###################################################################################################

waLBerla_compile_test( NAME PeriodicIntersect FILES PeriodicIntersect.cpp DEPENDS core blockforest )
waLBerla_execute_test( NAME PeriodicIntersect )
