
waLBerla_link_files_to_builddir( "*.dat" )


add_executable ( NonUniformGridBenchmark NonUniformGrid.cpp )
target_link_libraries ( NonUniformGridBenchmark walberla::blockforest walberla::boundary walberla::core
        walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop
        walberla::vtk walberla::sqlite )