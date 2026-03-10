if( WALBERLA_ENABLE_SWEEPGEN )
    add_subdirectory( DoubleShearLayer )
    add_subdirectory( FlowAroundSphere )
    add_subdirectory( ParallelPlates )
endif()

if( WALBERLA_BUILD_WITH_CODEGEN )
    add_subdirectory( MeshRefinementExample )
endif()

add_subdirectory( TimingPool )
