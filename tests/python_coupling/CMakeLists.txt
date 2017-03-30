###################################################################################################
#
# Tests for timeloop module
#
###################################################################################################

waLBerla_link_files_to_builddir( *.py )

if (WALBERLA_BUILD_WITH_PYTHON)
    waLBerla_compile_test( FILES ConfigFromPythonTest.cpp DEPENDS blockforest field )
    waLBerla_execute_test( NAME ConfigFromPythonTest 
                           COMMAND $<TARGET_FILE:ConfigFromPythonTest> ${CMAKE_CURRENT_SOURCE_DIR}/ConfigFromPythonTest.py )
                           
    waLBerla_compile_test( FILES CallbackTest.cpp  DEPENDS blockforest field )
    waLBerla_execute_test( NAME CallbackTest
                          COMMAND $<TARGET_FILE:CallbackTest> ${CMAKE_CURRENT_SOURCE_DIR}/CallbackTest.py ) 

    waLBerla_compile_test( FILES FieldExportTest.cpp  DEPENDS blockforest field )
    waLBerla_execute_test( NAME FieldExportTest
                           COMMAND $<TARGET_FILE:FieldExportTest> ${CMAKE_CURRENT_SOURCE_DIR}/FieldExportTest.py )
endif()
