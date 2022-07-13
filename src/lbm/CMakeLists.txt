###################################################################################################
#
# Module lbm
#
###################################################################################################

add_library( lbm )
target_link_libraries( lbm
      PUBLIC
      blockforest
      boundary
      communication
      core
      domain_decomposition
      field
      geometry
      gui
      stencil
      timeloop
      vtk
      )
target_sources( lbm
      PRIVATE
      PerformanceEvaluation.h
      BlockForestEvaluation.h
      IntelCompilerOptimization.h
      PerformanceLogger.h
      MassEvaluation.h
      )

add_subdirectory( vtk )
add_subdirectory( blockforest )
add_subdirectory( sweeps )
add_subdirectory( communication )
add_subdirectory( field )
add_subdirectory( free_surface )
add_subdirectory( refinement )
add_subdirectory( gui )
add_subdirectory( boundary )
add_subdirectory( mrt )
add_subdirectory( evaluations )
add_subdirectory( srt )
add_subdirectory( geometry )
add_subdirectory( trt )
add_subdirectory( cumulant )
add_subdirectory( lattice_model )
add_subdirectory( inplace_streaming )

###################################################################################################

