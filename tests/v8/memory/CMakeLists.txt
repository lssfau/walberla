
set( TestApp TestBufferSystem )

waLBerla_add_test_executable( ${TestApp} ${TestApp}.cpp )
target_link_libraries( ${TestApp} walberla::v8 )
walberla_v8_add_tests(
    ${TestApp}
    IDS
    testBufferIndexing
    testBufferIndexingWithPadding
    testConstCorrectness
)

set( TestApp TestAllocators )

waLBerla_add_test_executable( ${TestApp} ${TestApp}.cpp )

if( $CACHE{WALBERLA_BUILD_WITH_GPU_SUPPORT} )
    walberla_set_gpu_language( ${TestApp}.cu )
    target_sources( ${TestApp} PRIVATE ${TestApp}.cu )
endif()

target_link_libraries( ${TestApp} walberla::v8 )

walberla_v8_add_tests(
    ${TestApp}
    IDS
    testAlignment
    testAlignmentMT
    testOffsetAlignment
    testInterchange
    testBadAlloc
    testDefaultAndCopy
    testRebindAndConvert
    testSTLVector
    testSTLMap
)

set( TestApp TestField )

waLBerla_add_test_executable( ${TestApp} ${TestApp}.cpp )
target_link_libraries( ${TestApp} walberla::v8 )
walberla_v8_add_tests(
    ${TestApp}
    IDS
    testFieldCreation
    testConstCorrectness
    testFieldView
    testScalarFieldView
    testHostLineAlignment
)

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_HIP )

    set( TestApp TestFieldGPU )
    
    walberla_set_gpu_language( ${TestApp}.cu )
    waLBerla_add_test_executable( ${TestApp} ${TestApp}.cu )
    target_link_libraries( ${TestApp} walberla::v8 )
    walberla_v8_add_tests(
        ${TestApp}
        IDS
        testConstCorrectness
        testFieldViewUnified
        testCopyDeviceToPinned
        testScalarFieldView
        testDeviceLineAlignment
    )
endif()
