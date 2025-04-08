###################################################################################################
#
# Module lbm
#
###################################################################################################

add_library( walberla_lbm )
add_library( walberla::lbm ALIAS walberla_lbm )
target_link_libraries( walberla_lbm
      PUBLIC
      walberla::blockforest
      walberla::boundary
      walberla::communication
      walberla::core
      walberla::domain_decomposition
      walberla::field
      walberla::geometry
      walberla::stencil
      walberla::timeloop
      walberla::vtk
      )
target_sources( walberla_lbm
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

