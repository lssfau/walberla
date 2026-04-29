if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_HIP )
    set( TestApp TestGPUSweeper )

    walberla_set_gpu_language( ${TestApp}.cu )
    waLBerla_add_test_executable( ${TestApp} ${TestApp}.cu )
    target_link_libraries( ${TestApp} walberla::v8 )
    walberla_v8_add_tests(
        ${TestApp}
        IDS
        testForAllCells
        testForAllCellsStreams
        testForAllCellsStreamsConcurrency
    )
endif()