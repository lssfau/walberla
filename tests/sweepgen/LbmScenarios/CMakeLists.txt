
set( _testLanguages CXX )

if( $CACHE{WALBERLA_BUILD_WITH_CUDA} )
    list(APPEND _testLanguages CUDA )
endif()

if( $CACHE{WALBERLA_BUILD_WITH_HIP} )
    list(APPEND _testLanguages HIP )
endif()

set( TestIDs 
    FullyPeriodic
    MirroredHalfChannel
    FreeSlipPipe
)

set( _testExecutableTemplate TestLbmScenarios.in.cpp )

foreach(_testLanguage ${_testLanguages})
    
    set( _testExecutable ${CMAKE_CURRENT_BINARY_DIR}/TestLbmScenarios.${_testLanguage}.cpp )
    set( _testTarget TestLbmScenarios_${_testLanguage} )

    configure_file( ${_testExecutableTemplate} ${_testExecutable} )
    set_source_files_properties( ${_testExecutable} PROPERTIES LANGUAGE ${_testLanguage} )
    waLBerla_add_test_executable( ${_testTarget} ${_testExecutable} )
    target_sources(
        ${_testTarget} PRIVATE
        FILE_SET HEADERS
        FILES FullyPeriodic.hpp
              MirroredHalfChannel.hpp
    )
    
    walberla_generate_sources( ${_testTarget} 
        SCRIPTS LbmAlgorithms.py
        LANGUAGE ${_testLanguage}
        SCRIPT_ARGS --language=${_testLanguage} 
    )
    
    target_link_libraries( ${_testTarget} PRIVATE walberla::v8 )

    walberla_v8_add_tests(
        ${_testTarget}
        IDS ${TestIDs}
        PROCESSES 1 2 4
    )    
endforeach()
