add_executable( Test1D2DFields Test1D2DFields.cpp )
walberla_generate_sources( Test1D2DFields SCRIPTS Kernels1D2D.py )
target_link_libraries( Test1D2DFields PRIVATE walberla::core walberla::blockforest walberla::field walberla::experimental )
add_dependencies( waLBerlaTestsuite Test1D2DFields )

waLBerla_execute_test( NAME Test1D2DFields COMMAND Test1D2DFields )

add_executable( TestBlockforestGeometry TestBlockforestGeometry.cpp )
walberla_generate_sources( TestBlockforestGeometry SCRIPTS GeometryKernels.py )
target_link_libraries( TestBlockforestGeometry PRIVATE walberla::core walberla::blockforest walberla::field walberla::experimental )

add_dependencies( waLBerlaTestsuite TestBlockforestGeometry )

waLBerla_execute_test( NAME TestBlockforestGeometry COMMAND TestBlockforestGeometry )


add_subdirectory( BasicLbmScenarios )
add_subdirectory( boundaries )
