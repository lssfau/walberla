
# List modules here that should not be built
set ( excludedModules  )



# Generate the list of modulesToBuild ( all subfolders minus excludedModules )
file (GLOB allFilesAndDirs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} * )
list_minus ( modulesToBuild LIST1 ${allFilesAndDirs} LIST2 ${excludedModules} "CMakeFiles" )

# Activate modules
foreach ( mod ${modulesToBuild} )  
    if( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${mod} )       
       add_subdirectory( ${mod} ) 
    endif()
endforeach()       


# Write out a module statistic file that accumulates information about all 
# modules that have been added via waLBerla_add_module
waLBerla_write_module_statistics ( "${walberla_BINARY_DIR}/moduleStatistics.json" )


