###################################################################################################
#
# Tests for walberla::boundary   module
#
###################################################################################################


waLBerla_add_test_executable ( BoundaryHandling BoundaryHandling.cpp )
target_link_libraries ( BoundaryHandling walberla::core walberla::boundary )
waLBerla_execute_test( NAME BoundaryHandling )

if (WALBERLA_BUILD_WITH_PYTHON)

    waLBerla_link_files_to_builddir( *.py )

    waLBerla_add_test_executable ( TestShiftedPeriodicity TestShiftedPeriodicity.cpp )
    target_link_libraries ( TestShiftedPeriodicity walberla::boundary walberla::core walberla::blockforest
            walberla::field walberla::python_coupling )
    waLBerla_execute_test( NAME TestShiftedPeriodicity1 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py )
    waLBerla_execute_test( NAME TestShiftedPeriodicity2 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py PROCESSES 2 )
    waLBerla_execute_test( NAME TestShiftedPeriodicity4 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py PROCESSES 4 )

endif()