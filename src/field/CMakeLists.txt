
###################################################################################################
#
# Module field ( lattice data structure )
#
###################################################################################################

OPTION ( WALBERLA_FIELD_MONITORED_ACCESS "Enable (slow) monitoring mechanism for Fields" )
OPTION ( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION "Enables/Disables thread-safe field allocation (do NOT disable if you are not absolutely sure what you are doing!)" ON )
mark_as_advanced( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION )

configure_file ( CMakeDefs.in.h  CMakeDefs.h )

add_library( field )
target_link_libraries( field PUBLIC blockforest core communication domain_decomposition stencil vtk )
target_sources( field
      PRIVATE
      SwapableCompare.h
      AccuracyEvaluation.h
      FlagFunctions.h
      AccuracyEvaluationLinePlot.h
      VolumetricFlowRateEvaluation.h
      GhostLayerField.h
      GhostRegions.h
      AddToStorage.h
      EvaluationFilter.h
      FlagField.h
      Gather.h
      Layout.h
      SymmetryCheck.h
      FileIO.impl.h
      Field.impl.h
      Field.h
      FileIO.h
      FlagField.impl.h
      CMakeDefs.in.h
      CellCounter.h
      Printers.impl.h
      MassEvaluation.h
      GhostLayerField.impl.h
      Printers.h
      FieldClone.h
      StabilityChecker.h
      FlagUID.h
      )

add_subdirectory( adaptors )
add_subdirectory( allocation )
add_subdirectory( blockforest )
add_subdirectory( communication )
add_subdirectory( distributors )
add_subdirectory( interpolators )
add_subdirectory( iterators )
add_subdirectory( refinement )
add_subdirectory( vtk )