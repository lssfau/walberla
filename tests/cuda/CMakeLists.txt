###################################################################################################
#
# Tests for cuda
#
###################################################################################################

waLBerla_compile_test( FILES communication/GPUPackInfoTest.cpp DEPENDS blockforest )
waLBerla_execute_test( NAME  GPUPackInfoTest )

waLBerla_compile_test( FILES communication/GPUPackInfoCommunicationTest.cpp DEPENDS domain_decomposition blockforest stencil )
waLBerla_execute_test( NAME  GPUPackInfoCommunicationTest )

waLBerla_compile_test( FILES FieldTransferTest )
waLBerla_execute_test( NAME  FieldTransferTest )

waLBerla_compile_test( FILES SimpleKernelTest.cpp Kernels.cu DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME  SimpleKernelTest )

waLBerla_compile_test( FILES FieldIndexing3DTest.cpp FieldIndexing3DTest.cu )
waLBerla_execute_test( NAME  FieldIndexing3DTest )


waLBerla_python_file_generates(codegen/CudaJacobiKernel.py
        CudaJacobiKernel2D.cu CudaJacobiKernel2D.h
        CudaJacobiKernel3D.cu CudaJacobiKernel3D.h)
waLBerla_compile_test( FILES codegen/CodegenJacobiGPU.cpp
                             codegen/CudaJacobiKernel.py
                       DEPENDS blockforest timeloop gui )
waLBerla_execute_test( NAME CodegenJacobiGPU )


# The following tests work only for CUDA enabled MPI
waLBerla_compile_test( FILES communication/CommTest )
#waLBerla_execute_test( NAME  CommTest PROCESSES 2)

waLBerla_compile_test( FILES CudaMPI DEPENDS blockforest timeloop gui )
#waLBerla_execute_test( NAME  CudaMPI )

waLBerla_compile_test( FILES AlignmentTest.cpp DEPENDS blockforest timeloop )

waLBerla_python_file_generates(codegen/MicroBenchmarkGpuLbm.py MicroBenchmarkStreamKernel.cu MicroBenchmarkCopyKernel.cu)
waLBerla_compile_test( FILES codegen/MicroBenchmarkGpuLbm.cpp codegen/MicroBenchmarkGpuLbm.py)
