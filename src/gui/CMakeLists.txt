
# ---------------------------- GUI Module -----------------------------------------
#
# Directory Structure
#    /   Module base directory contains all files that can always be compiled
#        regardless if Qt was found - Qt dependend regions are marked with 
#        ifdef WALBERLA_ENABLE_GUI which is set from CMakeLists.txt
#   /BlockSliceView
#   /BlockView3D
#   /MainWindow
#         these 3 folders contain own code
#
#   /extern/Qt3D
#         source of Qt3D library, is needed until Qt 5 is released, which will 
#         contain the 3D library; until then, the needed parts of Qt3D are compiled 
#         as part of the the gui module
#   /extern/QtGradientEditor
#         Widget taken from Qt Designer - is needed for editing colormaps
#----------------------------------------------------------------------------------


# Sources from the gui base directory are always added, the parts including Qt Libs
# are marked with ifdefs
file( GLOB gui_SOURCES   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp   )
file( GLOB header        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h     )
file( GLOB cmake_rules   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cmake )
list( APPEND gui_SOURCES ${header} )
list( APPEND gui_SOURCES ${cmake_rules} )

if( WIN32 )
    SET   ( WALBERLA_GUI_USE_3D_BLOCKVIEW OFF CACHE INTERNAL "GUI: Use 3D BlockView - not supported on Windows!" FORCE )
else()
    OPTION( WALBERLA_GUI_USE_3D_BLOCKVIEW "GUI: Use 3D BlockView - not supported on Windows!" ON )
endif()


configure_file ( CMakeDefs.in.h  CMakeDefs.h )



# All other source files are only added if gui was activated
if( WALBERLA_ENABLE_GUI )
    include( ${CMAKE_CURRENT_SOURCE_DIR}/AddFilesQt.cmake ) 

    # -----------  Add Subfolders   --------------------------------------------
    
    include_directories (${CMAKE_CURRENT_BINARY_DIR} )
    set ( SUBFOLDERS ${CMAKE_CURRENT_SOURCE_DIR}/BlockSliceView 
                     ${CMAKE_CURRENT_SOURCE_DIR}/BlockDataDisplay
                     ${CMAKE_CURRENT_SOURCE_DIR}/BlockViewText
                     ${CMAKE_CURRENT_SOURCE_DIR}/MainWindow )

    if ( WALBERLA_GUI_USE_3D_BLOCKVIEW )
        list( APPEND SUBFOLDERS ${CMAKE_CURRENT_SOURCE_DIR}/BlockView3D )
    endif() 

    include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${SUBFOLDERS}  )
    add_files_qt(gui_SOURCES ${SUBFOLDERS} ) 
    
    
    # --------------- Extern GUI Packages---------------------------------------
    
    # Qt3D
    if ( WALBERLA_GUI_USE_3D_BLOCKVIEW )
        list_subdirectories (QT3D_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/extern/Qt3D )
        include_directories( SYSTEM ${QT3D_FOLDER} )
        add_files_qt(externSourceFiles ${QT3D_FOLDER} SOURCE_GROUP "extern\\\\Qt3D" )
    endif()

        
    # QtGradientEditor
    include_directories( SYSTEM ${QT3D_INCLUDE} 
                         ${CMAKE_CURRENT_SOURCE_DIR}/extern/QtGradientEditor )
    
    add_files_qt(externSourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/extern/QtGradientEditor SOURCE_GROUP "extern\\\\QtGradientEditor" )
    
    # Prevent warnings from extern files
    if ( WALBERLA_CXX_COMPILER_IS_MSVC )
        set_source_files_properties( ${externSourceFiles} PROPERTIES COMPILE_FLAGS "-W0" )
    else()
        set_source_files_properties( ${externSourceFiles} PROPERTIES COMPILE_FLAGS "-w" )
    endif()
    # add the extern files to module sources
    list(APPEND gui_SOURCES ${externSourceFiles})
    
endif(WALBERLA_ENABLE_GUI)


add_library( gui )
target_link_libraries( gui PUBLIC core domain_decomposition blockforest field timeloop stencil )
target_sources( gui
      PRIVATE
      PropertyTree.h
      PropertyTree.impl.h
      PropertyTree.cpp
      Gui.h
      CMakeDefs.in.h
      Gui.cpp
      )

#TODO: this needs some additional refactoring. Does not resamble the mechanism from above!
if( WALBERLA_ENABLE_GUI )
    add_subdirectory( MainWindow )
    add_subdirectory( BlockViewText )
    add_subdirectory( BlockView3D )
    add_subdirectory( BlockSliceView )
endif()
