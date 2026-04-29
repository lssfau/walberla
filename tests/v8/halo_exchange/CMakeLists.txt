
set( TestApp TestHaloExchange )

waLBerla_add_test_executable( ${TestApp} ${TestApp}.cpp )
target_link_libraries( ${TestApp} walberla::v8 )
walberla_v8_add_tests(
    ${TestApp}
    IDS
    testHaloExchangeInterfacesHost
    testHostFieldPackInfo
    testStreamPullPackInfo
    PROCESSES 1 2 4
)

if( $CACHE{WALBERLA_BUILD_WITH_GPU_SUPPORT} )
    walberla_set_gpu_language( ${TestApp}Device.cu )
    target_sources( ${TestApp} PRIVATE ${TestApp}Device.cu )
    
    walberla_v8_add_tests(
        ${TestApp}
        IDS
        testHaloExchangeInterfacesDevice
        testDeviceFieldPackInfo
        testDeviceStreamPullPackInfo
        PROCESSES 1 2 4
    )
endif()
