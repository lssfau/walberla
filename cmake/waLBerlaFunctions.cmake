include ( waLBerlaHelperFunctions         )
include ( waLBerlaModuleDependencySystem  )


set ( WALBERLA_GLOB_FILES *.cpp
                          *.c
                          *.h
                          *.cu
      CACHE INTERNAL "File endings to glob for source files" )


#######################################################################################################################
#
# Compiles an application either from given source files, or otherwise globs all files in the current folder.
# The application is linked against all waLBerla modules that are listed after DEPENDS
#
#
#   NAME    [required]    Name of application
#   GROUP   [optional]    IDE group name (e.g. VS)
#   DEPENDS [required]    list of modules, that this application depends on
#   FILES   [optional]    list of all source and header files belonging to that application
#                         if this is not given, all source and header files in the directory are added.
#                         Careful: when file was added or deleted, cmake has to be run again
#
#  Example:
#     waLBerla_compile_app ( NAME myApp DEPENDS core field lbm/boundary )
#######################################################################################################################

function ( waLBerla_add_executable )
    set( options )
    set( oneValueArgs NAME GROUP)
    set( multiValueArgs FILES DEPENDS)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( NOT ARG_NAME )
      message ( FATAL_ERROR "waLBerla_add_executable called without a NAME" )
    endif()

    # Skip this app, if it depends on modules that have not been built ( because they for example depend on PE)
    foreach ( depMod ${ARG_DEPENDS} )
        get_module_library_name ( depModLibraryName ${depMod} )
        if( NOT TARGET ${depModLibraryName} )
            if( WALBERLA_DEPS_ERROR )
               message( FATAL_ERROR "Module ${depMod} is missing to build target ${ARG_NAME}" )
            endif()
            if( WALBERLA_LOG_SKIPPED )
               message ( STATUS "Skipping ${ARG_NAME} since dependent module ${depMod} was not built" )
            endif()
            return()
        endif()
    endforeach()

    # Take source files either from parameter or search all source files
    set( sourceFiles )
    if ( ARG_FILES )
        foreach( sourceFile ${ARG_FILES} )
           get_filename_component( sourceFile ${sourceFile} ABSOLUTE )
           list( APPEND sourceFiles ${sourceFile} )
        endforeach()
        if ( WALBERLA_GROUP_FILES )
           group_files( "Source Files" FILES ${sourceFiles} )
        endif ( )
    else()
        file ( GLOB_RECURSE sourceFiles ${WALBERLA_GLOB_FILES} )  # Find all source files
        if ( WALBERLA_GROUP_FILES )
           file ( GLOB_RECURSE allFiles  "*" )  # Find all files
           list_minus ( otherFiles LIST1 ${allFiles} LIST2 ${sourceFiles} )
           group_files( "Source Files" FILES ${sourceFiles} )
           group_files( "Other Files" FILES ${otherFiles} )
        endif ( )
    endif()

    add_executable( ${ARG_NAME} ${sourceFiles} )
    target_link_libraries( ${ARG_NAME} ${WALBERLA_LINK_LIBRARIES_KEYWORD} ${ARG_DEPENDS} )
    target_link_libraries( ${ARG_NAME} ${WALBERLA_LINK_LIBRARIES_KEYWORD} ${SERVICE_LIBS} )
    set_property( TARGET ${ARG_NAME} PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD} )

    if( WALBERLA_GROUP_PROJECTS )
        if( NOT ARG_GROUP )
            set( ARG_GROUP "APPS" )
        endif()
        set_property( TARGET  ${ARG_NAME}  PROPERTY  FOLDER  ${ARG_GROUP} )
    endif()

endfunction ( waLBerla_add_executable )

#######################################################################################################################

#######################################################################################################################
#
# Adds a  waLBerla module test executable.
#
#######################################################################################################################

function ( waLBerla_compile_test )
    set( options )
    set( oneValueArgs NAME )
    set( multiValueArgs FILES DEPENDS )
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Module name is the directory relative to WALBERLA_MODULE_DIRS
    get_current_module_name ( )

    # Filename of first source file is used as name for testcase if no name was given
    if( NOT ARG_NAME )
        list( GET ARG_FILES 0 ARG_NAME )
        get_filename_component( ARG_NAME ${ARG_NAME} NAME_WE )
    endif()

    waLBerla_add_executable ( NAME ${ARG_NAME} GROUP "TESTS/${moduleName}"
                              DEPENDS ${ARG_DEPENDS} ${moduleName} FILES ${ARG_FILES} )


endfunction ( waLBerla_compile_test )
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

function ( waLBerla_link_files_to_builddir globExpression )

    # don't need links for in-source builds
    if( CMAKE_CURRENT_SOURCE_DIR STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
        return()
    endif()

    file( GLOB filesToLink RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${globExpression} )

    foreach( f ${filesToLink} )
        if( CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows" )
            configure_file( ${f} ${f} COPYONLY )
        else()
            execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${CMAKE_CURRENT_SOURCE_DIR}/${f}
                             ${CMAKE_CURRENT_BINARY_DIR}/${f} )
        endif()

    endforeach()

endfunction ( waLBerla_link_files_to_builddir )




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

function ( waLBerla_link_geometry_to_builddir globExpression )

   set(GEOMETRY_DIR ${walberla_SOURCE_DIR}/geometry)
   
   file(GLOB filesToLink "${GEOMETRY_DIR}/${globExpression}")
   
   # Create symlink link to the build directory (directly in ${CMAKE_CURRENT_BINARY_DIR})
   foreach(f ${filesToLink})
      get_filename_component(OBJ_FILENAME ${f} NAME)

      set(OBJ_BUILD_PATH "${CMAKE_CURRENT_BINARY_DIR}/${OBJ_FILENAME}")
      
      if( CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows" )
         configure_file(${f} ${OBJ_BUILD_PATH} COPYONLY)
      else()
        execute_process(   COMMAND ${CMAKE_COMMAND} -E create_symlink
                           ${f} ${OBJ_BUILD_PATH}
                        )
      endif()
   endforeach()

endfunction ( waLBerla_link_geometry_to_builddir )




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

function ( waLBerla_execute_test )

   set( options NO_MODULE_LABEL )
   set( oneValueArgs NAME PROCESSES )
   set( multiValueArgs COMMAND LABELS CONFIGURATIONS DEPENDS_ON_TARGETS )
   cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

   if( NOT ARG_NAME )
      message ( FATAL_ERROR "waLBerla_execute_test called without a NAME" )
   endif()

   if( NOT ARG_COMMAND AND NOT TARGET ${ARG_NAME} )
      if( WALBERLA_LOG_SKIPPED )
         message ( STATUS "Skipping test ${ARG_NAME} since the corresponding target is not built" )
      endif()
      return()
   endif()

   foreach( dependency_target ${ARG_DEPENDS_ON_TARGETS} )
      if( NOT TARGET ${dependency_target} )
         if( WALBERLA_LOG_SKIPPED )
            message ( STATUS "Skipping test ${ARG_NAME} since the target ${dependency_target} is not built" )
         endif()
         return()
      endif()
   endforeach( dependency_target )

   if( NOT ARG_PROCESSES )
      set ( numProcesses 1 )
   else()
      set ( numProcesses ${ARG_PROCESSES} )
   endif()

   if ( NOT ARG_COMMAND )
       set ( ARG_COMMAND $<TARGET_FILE:${ARG_NAME}> )
   endif()

   if( WALBERLA_BUILD_WITH_MPI )
      if( CMAKE_VERSION VERSION_LESS 3.10.0 )
	set ( MPIEXEC_EXECUTABLE ${MPIEXEC} )
      endif()
      list( INSERT  ARG_COMMAND  0  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numProcesses} ${MPIEXEC_PREFLAGS} )
   elseif( numProcesses GREATER 1 )
      return()
   endif()

   if( ARG_CONFIGURATIONS )
      add_test( NAME ${ARG_NAME} ${ARG_UNPARSED_ARGUMENTS} COMMAND ${ARG_COMMAND} CONFIGURATIONS ${ARG_CONFIGURATIONS} )
   else()
      add_test( NAME ${ARG_NAME} ${ARG_UNPARSED_ARGUMENTS} COMMAND ${ARG_COMMAND} )
   endif()

   if( ARG_NO_MODULE_LABEL )
      set_tests_properties ( ${ARG_NAME} PROPERTIES LABELS "${ARG_LABELS}" )
   else()
      get_current_module_name ( )
      set_tests_properties ( ${ARG_NAME} PROPERTIES LABELS "${moduleName} ${ARG_LABELS}" )
   endif()

   set_tests_properties ( ${ARG_NAME} PROPERTIES PROCESSORS ${numProcesses} )

endfunction ( waLBerla_execute_test )
#######################################################################################################################
