include ( waLBerlaHelperFunctions         )
include ( waLBerlaModuleDependencySystem  )


set ( WALBERLA_GLOB_FILES *.cpp
                          *.c
                          *.h
                          *.cu
      CACHE INTERNAL "File endings to glob for source files" )


#######################################################################################################################
#
# Creates a walberla module library
#
#
# Keywords:
#   DEPENDS [required]   list of modules, that this module depends on
#   FILES   [optional]   list of all source and header files belonging to that module
#                        if this is not given, all source and header files in the directory are added.
#                        Careful: when file was added or deleted, cmake has to be run again
#   EXCLUDE [optional]   Files that should not be included in the module (but are located in module directory).
#                        This makes only sense if FILES was not specified, and all files have been added automatically.
#   BUILD_ONLY_IF_FOUND  Before building the module test if all libraries specified here are availbable.
#   [optional]           This is done using the ${arg}_FOUND variable.
#                        Example: waLBerla_add_module( DEPENDS someModule BUILD_ONLY_IF_FOUND pe)
#                                 This module is only built if PE_FOUND is true.
#   OPTIONAL_DEPENDS     Lists modules, that this module might depend on. For example a module could depend on mesh_common if OpenMesh is
#   [optional]           available.
#
#######################################################################################################################

function ( waLBerla_add_module )
    set( options )
    set( oneValueArgs )
    set( multiValueArgs DEPENDS EXCLUDE FILES BUILD_ONLY_IF_FOUND OPTIONAL_DEPENDS )
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set( ALL_DEPENDENCIES ${ARG_DEPENDS} ${ARG_OPTIONAL_DEPENDS})
    # Module name is the directory relative to WALBERLA_MODULE_DIRS
    get_current_module_name ( moduleName )
    get_module_library_name ( moduleLibraryName ${moduleName} )

    # Test if all required libraries are available
    # this is detected using the _FOUND variable
    foreach ( externalName ${ARG_BUILD_ONLY_IF_FOUND} )
        string( TOUPPER ${externalName} externalNameUpper )
        if ( NOT ${externalNameUpper}_FOUND )
            message ( STATUS "Module ${moduleName} is not built, because ${externalName} not available" )
            return()
        endif()
    endforeach()

    # Take source files either from parameter or search all source files
    file ( GLOB_RECURSE allFiles "*" )  # Find all files
    if ( ARG_FILES )
        foreach( sourceFile ${ARG_FILES} )
           get_filename_component( sourceFile ${sourceFile} ABSOLUTE )
           list( APPEND sourceFiles ${sourceFile} )
        endforeach()
    else()
        file ( GLOB_RECURSE sourceFiles ${WALBERLA_GLOB_FILES} )  # Find all source files
    endif()

    # Remove exclude files from the sources
    if ( ARG_EXCLUDE )
        foreach ( fileToExclude ${ARG_EXCLUDE} )
            get_filename_component( fileToExclude ${fileToExclude} ABSOLUTE )
            list (REMOVE_ITEM sourceFiles ${fileToExclude} )
        endforeach()
    endif()

    list_minus ( otherFiles LIST1 ${allFiles} LIST2 ${sourceFiles} )
    set_source_files_properties( ${otherFiles} PROPERTIES HEADER_FILE_ONLY ON )

    if ( WALBERLA_GROUP_FILES )
       group_files( "Other Files"  FILES ${otherFiles}  )
       group_files( "Source Files" FILES ${sourceFiles} )
    endif ( )

    # Dependency Check
    check_dependencies( missingDeps additionalDeps FILES ${sourceFiles} EXPECTED_DEPS ${ALL_DEPENDENCIES} ${moduleName} )
    if ( missingDeps )
        message ( WARNING "The module ${moduleName} depends on ${missingDeps} which are not listed as dependencies!" )
    endif()

    set( hasSourceFiles FALSE )
 	foreach ( sourceFile ${sourceFiles} )
 	   if ( ${sourceFile} MATCHES "\\.(c|cpp|cu)" )
 	      set( hasSourceFiles TRUE )
 	   endif( )
 	endforeach( )

    if ( hasSourceFiles )
       add_library( ${moduleLibraryName} STATIC ${sourceFiles} ${otherFiles} )
 	else( )
       add_custom_target( ${moduleLibraryName} SOURCES ${sourceFiles} ${otherFiles} )  # dummy IDE target
 	endif( )

    waLBerla_register_dependency ( ${moduleName} ${ARG_DEPENDS} )

    set_property( TARGET ${moduleName} PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD} )

    # This property is needed for visual studio to group modules together
    if( WALBERLA_GROUP_PROJECTS )
       set_property( TARGET  ${moduleLibraryName}  PROPERTY  FOLDER  "SRC" )
    endif()

    # Install rule for library
    get_target_property( module_type ${moduleLibraryName} TYPE )
    if( ${module_type} MATCHES LIBRARY )
       install ( TARGETS ${moduleLibraryName}  RUNTIME DESTINATION bin
                                               LIBRARY DESTINATION lib
                                               ARCHIVE DESTINATION lib )
    endif( )

    # Install rule for header
    install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
              DESTINATION "walberla/${moduleName}"
              FILES_MATCHING PATTERN "*.h"
              PATTERN "*.in.h"     EXCLUDE
              PATTERN "CMakeFiles" EXCLUDE )

    # Install generated headers if necessary
    if ( NOT ${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
        install ( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
                  DESTINATION "walberla/${moduleName}"
                  FILES_MATCHING PATTERN "*.h"
                  PATTERN "*.in.h"     EXCLUDE
                  PATTERN "CMakeFiles" EXCLUDE )
    endif()


    # Report module statistics - which is later on written out to a json file
    # This file is used in the doxygen documentation to display a module graph
    #waLBerla_module_statistics ( FILES ${sourceFiles} DEPENDS ${ARG_DEPENDS} )

endfunction ( waLBerla_add_module )
#######################################################################################################################





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

    target_link_modules  ( ${ARG_NAME} ${ARG_DEPENDS}  )
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
    get_current_module_name ( moduleName )

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
      get_current_module_name ( moduleName  )
      set_tests_properties ( ${ARG_NAME} PROPERTIES LABELS "${moduleName} ${ARG_LABELS}" )
   endif()

   set_tests_properties ( ${ARG_NAME} PROPERTIES PROCESSORS ${numProcesses} )

endfunction ( waLBerla_execute_test )
#######################################################################################################################
