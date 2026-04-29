
set( TestApp TestSweepGeneration )

waLBerla_add_test_executable( ${TestApp} ${TestApp}.cpp )
walberla_generate_sources( ${TestApp} SCRIPTS ${TestApp}.py )
target_link_libraries( ${TestApp} PRIVATE walberla::core walberla::blockforest walberla::field walberla::v8 )

walberla_v8_add_tests(
    ${TestApp}
    IDS
    test1DField
    test2DField
    testV8Fields
    testV8FieldSwaps
    test1DFieldV8
    test2DFieldV8
)

add_executable( TestBlockforestGeometry TestBlockforestGeometry.cpp )
walberla_generate_sources( TestBlockforestGeometry SCRIPTS GeometryKernels.py )
target_link_libraries( TestBlockforestGeometry PRIVATE walberla::core walberla::blockforest walberla::field walberla::v8 )

add_dependencies( waLBerlaTestsuite TestBlockforestGeometry )

waLBerla_execute_test( NAME TestBlockforestGeometry COMMAND TestBlockforestGeometry )


add_subdirectory( BasicLbmScenarios )
add_subdirectory( boundaries )
