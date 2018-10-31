
waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME UniformGridBenchmarkGPU
                          FILES UniformGridGPU.cpp UniformGridGPU.gen.py
                          DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk )
