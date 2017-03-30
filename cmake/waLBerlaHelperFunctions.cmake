
#######################################################################################################################
#
# Function to add a flag to the given flag string
#
# param _VAR:  flag string
# param _FLAG: flag to add
# param ARGV0: optional boolean value if the flag should be added
#
#######################################################################################################################
FUNCTION ( add_flag  _VAR  _FLAG )
   if ( ARGC EQUAL 0 OR ARGV0 )
      set ( ${_VAR} "${${_VAR}} ${_FLAG}" PARENT_SCOPE )
   endif ( )
ENDFUNCTION ( add_flag )
#######################################################################################################################




#######################################################################################################################
#
# Subtracts list2 from list 2
# 
# Keywords:
#     LIST1   first list
#     LIST2   second list
# Example:
#     list_minus ( result LIST1 "entry1" "entry2" "entry3" LIST2 "entry1" "entry3" )
#      -> result has a single entry: "entry2"
#     
#######################################################################################################################
function( list_minus resultOut )
    set( options )
    set( oneValueArgs )
    set( multiValueArgs LIST1 LIST2 )
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
   
    set ( result ${ARG_LIST1} )
    foreach ( secondEntry ${ARG_LIST2} )
        list (REMOVE_ITEM result ${secondEntry} )
    endforeach()
    
    set ( ${resultOut} ${result} PARENT_SCOPE )
    
endfunction( list_minus )    
#######################################################################################################################




#######################################################################################################################
#
# Adds waLBerla to cmake package registry, and sets internal cache variables 
# ( include paths, compile definitions, compiler flags ... ) 
# 
# Which can be read by projects that add waLBerla as a subdirectory
# See for example waLBerla_import 
#    
#######################################################################################################################
function ( waLBerla_export )
    # Write the build directory to the cmake package config 
    # which is under  ~/.cmake/packages (Linux), or in the registry (Windows) 
    export( PACKAGE waLBerla )    
    
    # Copy a config file to the build directory which is executed when from another project
    # find_package( waLBerla ) is called 
    configure_file ( cmake/walberla-config-builddir.cmake
                     walberla-config.cmake @ONLY )
    
    # Export compiler flags
    set ( compilerList "CXX" "C" )
    foreach ( compiler ${compilerList} )
        set ( WALBERLA_INCLUDE_SYSTEM_FLAG_${compiler}  ${CMAKE_INCLUDE_SYSTEM_FLAG_${compiler}}
              CACHE INTERNAL "${compiler} include flags for walberla" )
        set ( WALBERLA_${compiler}_FLAGS ${CMAKE_${compiler}_FLAGS} 
              CACHE INTERNAL "${compiler} flags for walberla" )
              
        foreach ( buildType ${CMAKE_CONFIGURATION_TYPES} )        
            set ( WALBERLA_${compiler}_FLAGS_${buildType} ${CMAKE_${compiler}_FLAGS_${buildType}} 
                  CACHE INTERNAL "${compiler} flags for walberla" )
        endforeach()
    endforeach()    
    
    # Export linker flags
    set ( WALBERLA_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} CACHE INTERNAL "walberla linker flags")
    foreach ( buildType ${CMAKE_CONFIGURATION_TYPES} )        
        set ( WALBERLA_EXE_LINKER_FLAGS_${buildType} ${CMAKE_EXE_LINKER_FLAGS_${buildType}} 
                CACHE INTERNAL "walberla linker flags for ${buildType}" )
    endforeach()
    
    # Export service libs
    set ( WALBERLA_SERVICE_LIBS ${SERVICE_LIBS} CACHE INTERNAL "External Libraries necessary for waLBerla" )
    
    # Export compile definitions
    get_directory_property( WALBERLA_COMPILE_DEFINITIONS DIRECTORY ${walberla_SOURCE_DIR} COMPILE_DEFINITIONS )
    set ( WALBERLA_COMPILE_DEFINITIONS ${WALBERLA_COMPILE_DEFINITIONS} CACHE INTERNAL "waLBerla compile definitions" )

    # Export include paths
    get_directory_property( WALBERLA_INCLUDE_DIRS DIRECTORY ${walberla_SOURCE_DIR} INCLUDE_DIRECTORIES)
    set ( WALBERLA_INCLUDE_DIRS ${WALBERLA_INCLUDE_DIRS} CACHE INTERNAL "waLBerla include directories" )
    
    # Export link paths
    set ( WALBERLA_LINK_DIRS ${LINK_DIRS} CACHE INTERNAL "waLBerla link directories" )
endfunction( waLBerla_export)

#######################################################################################################################




#######################################################################################################################
#
# Reads variables that have been exported by waLBerla_export
# Useful when waLBerla was added as a subdirectory
#
#######################################################################################################################
function ( waLBerla_import )
    set( options NO_COMPILER_FLAGS NO_INCLUDE_DIRS NO_COMPILE_DEFINITIONS )
    set( oneValueArgs )
    set( multiValueArgs  )
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( NOT NO_INCLUDE_DIRS )
        foreach( directory ${WALBERLA_INCLUDE_DIRS} )
            if( ${directory} MATCHES "walberla/src$" )
                include_directories( ${directory} )
            else()
                include_directories( SYSTEM ${directory} )
            endif()
        endforeach()
    endif()
    
    # Import compile definitions
    if ( NOT NO_COMPILE_DEFINITIONS )
        set_property( DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS ${WALBERLA_COMPILE_DEFINITIONS} )
    endif()
        
    # Import compiler flags
    if ( NOT NO_COMPILER_FLAGS )
        set ( compilerList "CXX" "C" )
        foreach ( compiler ${compilerList} )
            set ( CMAKE_INCLUDE_SYSTEM_FLAG_${compiler}  ${WALBERLA_INCLUDE_SYSTEM_FLAG_${compiler}}  PARENT_SCOPE )
            set ( CMAKE_${compiler}_FLAGS                ${WALBERLA_${compiler}_FLAGS}                PARENT_SCOPE )
                  
            foreach ( buildType ${CMAKE_CONFIGURATION_TYPES} )        
                set ( CMAKE_${compiler}_FLAGS_${buildType} ${WALBERLA_${compiler}_FLAGS_${buildType}} PARENT_SCOPE )
            endforeach()
        endforeach()    
        
        # Import linker flags
        set ( CMAKE_EXE_LINKER_FLAGS ${WALBERLA_EXE_LINKER_FLAGS}  )
        foreach ( buildType ${CMAKE_CONFIGURATION_TYPES} )        
            set ( CMAKE_EXE_LINKER_FLAGS_${buildType} ${WALBERLA_EXE_LINKER_FLAGS_${buildType}} PARENT_SCOPE )
        endforeach() 
    endif()
    
    set( SERVICE_LIBS ${WALBERLA_SERVICE_LIBS} PARENT_SCOPE )
    
    link_directories( ${WALBERLA_LINK_DIRS} )
endfunction( waLBerla_import)
#######################################################################################################################


#######################################################################################################################
#
# Group Files for an IDE like VS
# 
# Keywords:
#     FILES     all files which are in a certain group
# Example:
#     file        ( GLOB_RECURSE sourceFiles "*.h" "*.c" )
#     group_files ( "Source Files" FILES ${sourceFiles} )
#     
#######################################################################################################################
function( group_files group )
    set( options )
    set( oneValueArgs)
    set( multiValueArgs FILES)
    cmake_parse_arguments( ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
    
    foreach( filename ${ARG_FILES} )
       file( RELATIVE_PATH rel ${CMAKE_CURRENT_SOURCE_DIR} ${filename})
       if( rel MATCHES ".*/.*" )
          # subdirectory found -> add to subdirectory
          get_filename_component(dir ${rel} PATH )
          # replace slashes by some underscores to use path as variable name 
          # -> you are allowed to use at max two underscores in your folder name
          string( REGEX REPLACE  "/" "___" subname ${dir} )
          list( APPEND subfiles${subname} ${filename} )
          list( APPEND subnames           ${subname}  )
       else()
          # no subdirectory found -> add to normal group
          list( APPEND files ${filename} )
       endif()
    endforeach()
    
    if( subnames )
       list( REMOVE_DUPLICATES subnames )
    endif()
    
    foreach( subname ${subnames} )
       string( REGEX REPLACE  "___" "\\\\" subgroup ${subname} )
       #message( STATUS " Group: ${group}\\${subgroup} and files: ${subfiles${subname}}" )
       source_group( "${group}\\${subgroup}" FILES ${subfiles${subname}} )
    endforeach()
    #message( STATUS " Group: ${group} and files: ${files}" )
    source_group( "${group}" FILES ${files} )
endfunction( group_files )
#######################################################################################################################



#######################################################################################################################
#
# Sets version number of waLberla
# 
# Example:
#     
#######################################################################################################################
function( set_version VERSION_MAJOR VERSION_PATCH )
  set( WALBERLA_MAJOR_VERSION ${VERSION_MAJOR} CACHE STRING "waLBerla major version" FORCE )
  set( WALBERLA_PATCH_LEVEL   ${VERSION_PATCH} CACHE STRING "waLBerla patch level" FORCE )
  set( WALBERLA_VERSION "${VERSION_MAJOR}.${VERSION_PATCH}" CACHE STRING "waLBerla version" FORCE )
  mark_as_advanced( WALBERLA_MAJOR_VERSION WALBERLA_PATCH_LEVEL )
endfunction()
#######################################################################################################################
