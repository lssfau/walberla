
set( TestIDs 
    FullyPeriodic
    MirroredHalfChannel
    FreeSlipPipe
)

add_executable( TestBasicLbmScenariosCPU TestBasicLbmScenarios.cpp )
walberla_generate_sources( TestBasicLbmScenariosCPU SCRIPTS LbmAlgorithms.py SCRIPT_ARGS --target=cpu )
target_link_libraries( TestBasicLbmScenariosCPU PRIVATE walberla::core walberla::blockforest walberla::field walberla::geometry walberla::experimental )
add_dependencies( waLBerlaTestsuite TestBasicLbmScenariosCPU )

foreach( TestID ${TestIDs} )
foreach( NumProcs 1 2 4 )
    waLBerla_execute_test(
        NAME "BasicLbmScenarios - CPU - ${NumProcs} Processes - ${TestID}"
        COMMAND $<TARGET_FILE:TestBasicLbmScenariosCPU> ${TestID} 
        PROCESSES ${NumProcs}
    )
endforeach()
endforeach()


if( $CACHE{WALBERLA_BUILD_WITH_CUDA} )
    set( _codegen_suffixes hpp cu )

    add_executable( TestBasicLbmScenariosCUDA TestBasicLbmScenarios.cpp )
    walberla_generate_sources( TestBasicLbmScenariosCUDA SCRIPTS LbmAlgorithms.py SCRIPT_ARGS --target=cuda FILE_EXTENSIONS ${_codegen_suffixes} )
    target_link_libraries( TestBasicLbmScenariosCUDA PRIVATE walberla::core walberla::blockforest walberla::field walberla::gpu walberla::geometry walberla::experimental CUDA::cudart )

    add_dependencies( waLBerlaTestsuite TestBasicLbmScenariosCUDA )

    foreach( TestID ${TestIDs} )
    foreach( NumProcs 1 2 4 )
        waLBerla_execute_test(
            NAME "BasicLbmScenarios - CUDA - ${NumProcs} Processes - ${TestID}"
            COMMAND $<TARGET_FILE:TestBasicLbmScenariosCUDA> ${TestID} 
            PROCESSES ${NumProcs}
        )
    endforeach()
    endforeach()
endif()

if( $CACHE{WALBERLA_BUILD_WITH_HIP} )
    set( _codegen_suffixes hpp hip )

    add_executable( TestBasicLbmScenariosHIP TestBasicLbmScenarios.cpp )
    walberla_generate_sources( TestBasicLbmScenariosHIP SCRIPTS LbmAlgorithms.py SCRIPT_ARGS --target=hip FILE_EXTENSIONS ${_codegen_suffixes} )
    target_link_libraries( TestBasicLbmScenariosHIP PRIVATE walberla::core walberla::blockforest walberla::field walberla::gpu walberla::geometry walberla::experimental )

    add_dependencies( waLBerlaTestsuite TestBasicLbmScenariosHIP )

    foreach( TestID ${TestIDs} )
        foreach( NumProcs 1 2 4 )
        waLBerla_execute_test(
            NAME "BasicLbmScenarios - HIP - ${NumProcs} Processes - ${TestID}"
            COMMAND $<TARGET_FILE:TestBasicLbmScenariosHIP> ${TestID} 
            PROCESSES ${NumProcs}
        )
    endforeach()
    endforeach()
endif()
