
waLBerla_link_files_to_builddir( "*.dat" )


waLBerla_add_executable ( NAME UniformGridBenchmark 
                          DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk sqlite )