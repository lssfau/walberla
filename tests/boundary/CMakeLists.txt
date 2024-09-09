###################################################################################################
#
# Tests for boundary module
#
###################################################################################################


waLBerla_compile_test( FILES BoundaryHandling.cpp )
waLBerla_execute_test( NAME BoundaryHandling )

if (WALBERLA_BUILD_WITH_PYTHON)

    waLBerla_link_files_to_builddir( *.py )

    waLBerla_compile_test( FILES TestShiftedPeriodicity.cpp DEPENDS blockforest field python_coupling )
    waLBerla_execute_test( NAME TestShiftedPeriodicity1 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py )
    waLBerla_execute_test( NAME TestShiftedPeriodicity2 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py PROCESSES 2 )
    waLBerla_execute_test( NAME TestShiftedPeriodicity4 COMMAND $<TARGET_FILE:TestShiftedPeriodicity> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetup.py PROCESSES 4 )

endif()