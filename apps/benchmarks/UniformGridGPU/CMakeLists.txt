
waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME UniformGridBenchmarkGPU
                          FILES UniformGridGPU.cpp UniformGridGPU_LatticeModel.cpp
                                UniformGridGPU_LbKernel.cu UniformGridGPU_NoSlip.cu UniformGridGPU_UBB.cu
                                UniformGridGPU_PackInfo.cu
                          DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk )
