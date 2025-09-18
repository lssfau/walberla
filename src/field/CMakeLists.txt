
###################################################################################################
#
# Module walberla::field ( lattice data structure )
#
###################################################################################################

OPTION ( WALBERLA_FIELD_MONITORED_ACCESS "Enable (slow) monitoring mechanism for Fields" )
OPTION( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION "Enables/Disables thread-safe walberla::field allocation (do NOT disable if you are not absolutely sure what you are doing!)" ON )
mark_as_advanced( WALBERLA_THREAD_SAFE_FIELD_ALLOCATION )

configure_file ( CMakeDefs.in.h  CMakeDefs.h )

add_library( walberla_field )
add_library( walberla::field ALIAS walberla_field )
target_link_libraries( walberla_field PUBLIC walberla::blockforest walberla::core walberla::communication walberla::domain_decomposition walberla::stencil walberla::vtk )
target_sources( walberla_field
        PUBLIC
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
        Field.h
        FileIO.h
        CMakeDefs.in.h
        CellCounter.h
        MassEvaluation.h
        Printers.h
        FieldClone.h
        StabilityChecker.h
        FlagUID.h
        FileIO.impl.h
        Field.impl.h
        FlagField.impl.h
        Printers.impl.h
        GhostLayerField.impl.h
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