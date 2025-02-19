
waLBerla_link_files_to_builddir( "*.dat" )


waLBerla_add_executable ( NAME UniformGridBenchmark
      DEPENDS walberla::blockforest walberla::boundary walberla::core walberla::domain_decomposition walberla::field walberla::lbm walberla::postprocessing walberla::timeloop walberla::vtk walberla::sqlite )