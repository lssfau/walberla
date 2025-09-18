###################################################################################################
#
# Tests for python coupling
#
###################################################################################################

waLBerla_link_files_to_builddir( *.py )

if (WALBERLA_BUILD_WITH_PYTHON)
    waLBerla_add_test_executable ( ConfigFromPythonTest ConfigFromPythonTest.cpp )
    target_link_libraries ( ConfigFromPythonTest walberla::blockforest walberla::field walberla::python_coupling )
    waLBerla_execute_test( NAME ConfigFromPythonTest 
                           COMMAND $<TARGET_FILE:ConfigFromPythonTest> ${CMAKE_CURRENT_SOURCE_DIR}/ConfigFromPythonTest.py )

    waLBerla_add_test_executable ( CallbackTest CallbackTest.cpp )
    target_link_libraries ( CallbackTest walberla::blockforest walberla::field walberla::python_coupling )
    waLBerla_execute_test( NAME CallbackTest
                          COMMAND $<TARGET_FILE:CallbackTest> ${CMAKE_CURRENT_SOURCE_DIR}/CallbackTest.py )

    waLBerla_add_test_executable ( FieldExportTest FieldExportTest.cpp )
    target_link_libraries ( FieldExportTest walberla::blockforest walberla::field walberla::python_coupling )
    waLBerla_execute_test( NAME FieldExportTest
                           COMMAND $<TARGET_FILE:FieldExportTest> ${CMAKE_CURRENT_SOURCE_DIR}/FieldExportTest.py )

    add_dependencies( waLBerlaTestsuite pythonModule )
    add_test( NAME PythonWalberlaTest
          COMMAND ${Python_EXECUTABLE} -m unittest discover waLBerla_tests/ WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python/)
    SET_TESTS_PROPERTIES(PythonWalberlaTest
          PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/apps/pythonmodule:$PYTHONPATH")
endif()
