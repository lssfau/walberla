
waLBerla_link_files_to_builddir( "*.dat" )


add_executable ( UniformGridBenchmark UniformGrid.cpp )
target_link_libraries ( UniformGridBenchmark walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk walberla::sqlite )