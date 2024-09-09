
# List modules here that should not be built
#set ( excludedModules core )
#
#
#
## Generate the list of modulesToBuild ( all subfolders minus excludedModules )
#file (GLOB allFilesAndDirs  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} * )
#list_minus ( modulesToBuild LIST1 ${allFilesAndDirs} LIST2 ${excludedModules} "CMakeFiles" )
#
## Activate modules
#foreach ( mod ${modulesToBuild} )
#    if( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${mod} )
#       add_subdirectory( ${mod} )
#    endif()
#endforeach()
#
#
## Write out a module statistic file that accumulates information about all
## modules that have been added via waLBerla_add_module
#waLBerla_write_module_statistics ( "${walberla_BINARY_DIR}/moduleStatistics.json" )

add_subdirectory( blockforest )
add_subdirectory( boundary )
add_subdirectory( communication )
add_subdirectory( core )
add_subdirectory( gpu )
add_subdirectory( domain_decomposition )
add_subdirectory( executiontree )
if ( WALBERLA_BUILD_WITH_FFT AND FFTW3_FOUND )
   add_subdirectory( fft )
endif()
add_subdirectory( field )
add_subdirectory( gather )
add_subdirectory( geometry )
add_subdirectory( gui )
add_subdirectory( lbm )
add_subdirectory( lbm_generated )
add_subdirectory( lbm_mesapd_coupling )
add_subdirectory( mesa_pd )
if( OPENMESH_FOUND )
   add_subdirectory( mesh )
   add_subdirectory( mesh_common )
endif()
add_subdirectory( pde )
add_subdirectory( pe )
add_subdirectory( pe_coupling )
add_subdirectory( postprocessing )
if (WALBERLA_BUILD_WITH_PYTHON)
   add_subdirectory(python_coupling)
endif ()
add_subdirectory( simd )
add_subdirectory( sqlite )
add_subdirectory( stencil )
add_subdirectory( timeloop )
add_subdirectory( vtk )

