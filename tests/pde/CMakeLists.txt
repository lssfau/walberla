#############################################################################################################################
#
# Tests for walberla::pde module
#
#############################################################################################################################


waLBerla_add_test_executable ( CGTest CGTest.cpp )
target_link_libraries ( CGTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME CGShortTest COMMAND $<TARGET_FILE:CGTest> --shortrun PROCESSES 8 )

waLBerla_add_test_executable ( JacobiTest JacobiTest.cpp )
target_link_libraries ( JacobiTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME JacobiShortTest COMMAND $<TARGET_FILE:JacobiTest> --shortrun PROCESSES 8 )

waLBerla_add_test_executable ( RBGSTest RBGSTest.cpp )
target_link_libraries ( RBGSTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME RBGSShortTest COMMAND $<TARGET_FILE:RBGSTest> --shortrun PROCESSES 8 )

waLBerla_add_test_executable ( SORTest SORTest.cpp )
target_link_libraries ( SORTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME SORShortTest COMMAND $<TARGET_FILE:SORTest> --shortrun PROCESSES 8 )

waLBerla_add_test_executable ( MGTest MGTest.cpp )
target_link_libraries ( MGTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME MGShortTest COMMAND $<TARGET_FILE:MGTest> --shortrun PROCESSES 8 )
waLBerla_execute_test( NAME MGTest COMMAND $<TARGET_FILE:MGTest> PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( MGConvergenceTest MGConvergenceTest.cpp )
target_link_libraries ( MGConvergenceTest walberla::blockforest walberla::timeloop walberla::vtk walberla::core walberla::pde )
waLBerla_execute_test( NAME MGConvergenceTest COMMAND $<TARGET_FILE:MGConvergenceTest> PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )

waLBerla_add_test_executable ( BoundaryTest BoundaryTest.cpp )
target_link_libraries ( BoundaryTest walberla::blockforest walberla::timeloop walberla::vtk walberla::boundary walberla::core walberla::pde )
waLBerla_execute_test( NAME BoundaryShortTest COMMAND $<TARGET_FILE:BoundaryTest> --shortrun PROCESSES 8 )
waLBerla_execute_test( NAME BoundaryTest COMMAND $<TARGET_FILE:BoundaryTest> PROCESSES 8 CONFIGURATIONS Release RelWithDbgInfo )
