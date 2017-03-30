#############################################################################################################################
#
# Tests for pde module
#
#############################################################################################################################


waLBerla_compile_test( FILES CGTest.cpp DEPENDS blockforest timeloop vtk )
waLBerla_execute_test( NAME CGShortTest COMMAND $<TARGET_FILE:CGTest> --shortrun PROCESSES 8 )

waLBerla_compile_test( FILES JacobiTest.cpp DEPENDS blockforest timeloop vtk )
waLBerla_execute_test( NAME JacobiShortTest COMMAND $<TARGET_FILE:JacobiTest> --shortrun PROCESSES 8 )

waLBerla_compile_test( FILES RBGSTest.cpp DEPENDS blockforest timeloop vtk )
waLBerla_execute_test( NAME RBGSShortTest COMMAND $<TARGET_FILE:RBGSTest> --shortrun PROCESSES 8 )

waLBerla_compile_test( FILES SORTest.cpp DEPENDS blockforest timeloop vtk )
waLBerla_execute_test( NAME SORShortTest COMMAND $<TARGET_FILE:SORTest> --shortrun PROCESSES 8 )

waLBerla_compile_test( FILES MGTest.cpp DEPENDS blockforest timeloop vtk )
waLBerla_execute_test( NAME MGShortTest COMMAND $<TARGET_FILE:MGTest> --shortrun PROCESSES 8 )
waLBerla_execute_test( NAME MGTest COMMAND $<TARGET_FILE:MGTest> PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )
