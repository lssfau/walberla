
waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_python_file_generates(UniformGridGPU.py
        UniformGridGPU_LatticeModel.cpp UniformGridGPU_LatticeModel.h
        UniformGridGPU_LbKernel.cu UniformGridGPU_LbKernel.h
        UniformGridGPU_NoSlip.cu UniformGridGPU_NoSlip.h
        UniformGridGPU_UBB.cu UniformGridGPU_UBB.h
        UniformGridGPU_PackInfo.cu UniformGridGPU_PackInfo.h
        )

waLBerla_add_executable ( NAME UniformGridBenchmarkGPU
                          FILES UniformGridGPU.cpp UniformGridGPU.py
                          DEPENDS blockforest boundary core cuda domain_decomposition field geometry timeloop vtk )
