###################################################################################################
#
# Tests for gather module
#
###################################################################################################


waLBerla_compile_test( FILES MPIGatherSchemeTest.cpp DEPENDS )
waLBerla_execute_test( NAME MPIGatherSchemeTest PROCESSES 7 )

waLBerla_compile_test( FILES CurveGatherTest.cpp DEPENDS )
waLBerla_execute_test( NAME CurveGatherTest PROCESSES 4 )

waLBerla_compile_test( FILES GatherSchemeTest.cpp DEPENDS )
waLBerla_execute_test( NAME GatherSchemeTest PROCESSES 4 )
