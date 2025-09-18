include(waLBerlaHelperFunctions)


set(WALBERLA_GLOB_FILES *.cpp
        *.c
        *.h
        *.cu
        CACHE INTERNAL "File endings to glob for source files")

#######################################################################################################################

#######################################################################################################################
#
# Adds a waLBerla test executable.
#
#######################################################################################################################

add_custom_target(waLBerlaTestsuite)

function(waLBerla_add_test_executable)
    add_executable(${ARGV})
    add_dependencies(waLBerlaTestsuite ${ARGV0})
endfunction()

#######################################################################################################################


#######################################################################################################################
#
# Links all files in current source dir matching a globbing expression to the build directory
#
# first parameter is glob expression
#
# Typical usage: link all parameter files in same folder as the binary was produced
#                Assuming the parameter files end with prm, write this to your CMakeLists in the
#                app or test folder:
#                waLBerla_link_files_to_builddir( "*.prm" )
#
# Symlinking works only under linux, on windows the files are copied. For in-source builds this does nothing.
#
#######################################################################################################################

function(waLBerla_link_files_to_builddir globExpression)

    # don't need links for in-source builds
    if (CMAKE_CURRENT_SOURCE_DIR STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
        return()
    endif ()

    file(GLOB filesToLink RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${globExpression})

    foreach (f ${filesToLink})
        if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
            configure_file(${f} ${f} COPYONLY)
        else ()
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${CMAKE_CURRENT_SOURCE_DIR}/${f}
                    ${CMAKE_CURRENT_BINARY_DIR}/${f})
        endif ()

    endforeach ()

endfunction(waLBerla_link_files_to_builddir)


#######################################################################################################################
#
# Links all files in PROJECT_ROOT_DIR/geometry matching a globbing expression to the build directory
#
# first parameter is glob expression
#
# Typical usage: Assuming the parameter files end with obj, write this to your CMakeLists in the
#                app or test folder:
#                waLBerla_link_geometry_to_builddir( "*.obj" )
#
# Symlinking works only under linux, on windows the files are copied.
#
#######################################################################################################################

function(waLBerla_link_geometry_to_builddir globExpression)

    set(GEOMETRY_DIR ${walberla_SOURCE_DIR}/geometry)

    file(GLOB filesToLink "${GEOMETRY_DIR}/${globExpression}")

    # Create symlink link to the build directory (directly in ${CMAKE_CURRENT_BINARY_DIR})
    foreach (f ${filesToLink})
        get_filename_component(OBJ_FILENAME ${f} NAME)

        set(OBJ_BUILD_PATH "${CMAKE_CURRENT_BINARY_DIR}/${OBJ_FILENAME}")

        if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
            configure_file(${f} ${OBJ_BUILD_PATH} COPYONLY)
        else ()
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    ${f} ${OBJ_BUILD_PATH}
            )
        endif ()
    endforeach ()

endfunction(waLBerla_link_geometry_to_builddir)


#######################################################################################################################
#
# Determine Module name using the current folder
#
# moduleFolder is the current source directory relative to a folder in WALBERLA_MODULE_DIRS
# The variable moduleName will be set in PARENT_SCOPE and is the first folder in WALBERLA_MODULE_DIRS
# Example:
#    If CMAKE_CURRENT_SOURCE_DIR is /src/core/field and /src/ is an element in WALBERLA_MODULE_DIRS,
#    then moduleName is "core"
#
#######################################################################################################################

function(waLBerla_current_module_name)
    foreach (moduleDir ${ARGN} ${WALBERLA_MODULE_DIRS})
        file(RELATIVE_PATH moduleFolder ${moduleDir} ${CMAKE_CURRENT_SOURCE_DIR})
        if (NOT ${moduleFolder} MATCHES "\\.\\./.*")
            #append / to make cmake_path also work with one directory only
            string(REGEX REPLACE "(.*)/.*" "\\1" moduleNameOut ${moduleFolder})
            set(moduleName walberla::${moduleNameOut} PARENT_SCOPE)
            return()
        endif ()
    endforeach ()

    message(WARNING "Called get_current_module_name, in a directory "
            "that is not a subdirectory of WALBERLA_MODULE_DIRS\n"
            "Module Dirs: ${WALBERLA_MODULE_DIRS} \n"
            "Current Dir: ${CMAKE_CURRENT_SOURCE_DIR}")


endfunction(waLBerla_current_module_name)


#######################################################################################################################
#
# Adds an executable to CTest.
#
# Wrapper around add_test, that handles test labels and parallel runs
# Adds the module name as test label, plus optional user defined labels.
#
#   NAME      [required] Name of test
#   PROCESSES [optional] Number of MPI processes, that are used to start this test.
#                        If walberla is built without MPI, and PROCESSES > 1, the test is not added
#                        Defaults to 1
#   COMMAND   [optional] The command that is executed. Use this to start with parameter files or other
#                        command line options.
#                        Defaults to $<TARGET_FILE:${NAME}> (for this syntax see cmake documentation of add_test )
#   LABELS    [optional] Additional test labels.
#
#######################################################################################################################

function(waLBerla_execute_test)

    set(options NO_MODULE_LABEL)
    set(oneValueArgs NAME PROCESSES)
    set(multiValueArgs COMMAND LABELS CONFIGURATIONS DEPENDS_ON_TARGETS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT ARG_NAME)
        message(FATAL_ERROR "waLBerla_execute_test called without a NAME")
    endif ()

    if (NOT ARG_COMMAND AND NOT TARGET ${ARG_NAME})
        if (WALBERLA_LOG_SKIPPED)
            message(STATUS "Skipping test ${ARG_NAME} since the corresponding target is not built")
        endif ()
        return()
    endif ()

    foreach (dependency_target ${ARG_DEPENDS_ON_TARGETS})
        if (NOT TARGET ${dependency_target})
            if (WALBERLA_LOG_SKIPPED)
                message(STATUS "Skipping test ${ARG_NAME} since the target ${dependency_target} is not built")
            endif ()
            return()
        endif ()
    endforeach (dependency_target)

    if (NOT ARG_PROCESSES)
        set(numProcesses 1)
    else ()
        set(numProcesses ${ARG_PROCESSES})
    endif ()

    if (NOT ARG_COMMAND)
        set(ARG_COMMAND $<TARGET_FILE:${ARG_NAME}>)
    endif ()

    if (WALBERLA_BUILD_WITH_MPI)
        if (CMAKE_VERSION VERSION_LESS 3.10.0)
            set(MPIEXEC_EXECUTABLE ${MPIEXEC})
        endif ()
        list(INSERT ARG_COMMAND 0 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numProcesses} ${MPIEXEC_PREFLAGS})
    elseif (numProcesses GREATER 1)
        return()
    endif ()

    if (ARG_CONFIGURATIONS)
        add_test(NAME ${ARG_NAME} ${ARG_UNPARSED_ARGUMENTS} COMMAND ${ARG_COMMAND} CONFIGURATIONS ${ARG_CONFIGURATIONS})
    else ()
        add_test(NAME ${ARG_NAME} ${ARG_UNPARSED_ARGUMENTS} COMMAND ${ARG_COMMAND})
    endif ()

    if (ARG_NO_MODULE_LABEL)
        set_tests_properties(${ARG_NAME} PROPERTIES LABELS "${ARG_LABELS}")
    else ()
        waLBerla_current_module_name ( )
        set_tests_properties(${ARG_NAME} PROPERTIES LABELS "${moduleName} ${ARG_LABELS}")
    endif ()

    set_tests_properties(${ARG_NAME} PROPERTIES PROCESSORS ${numProcesses})

endfunction(waLBerla_execute_test)
#######################################################################################################################
