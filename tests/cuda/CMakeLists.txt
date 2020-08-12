###################################################################################################
#
# Tests for cuda
#
###################################################################################################

if( WALBERLA_BUILD_WITH_CUDA )

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


waLBerla_generate_target_from_python(NAME CodegenJacobiGPUGeneratedCudaJacobiKernel FILE codegen/CudaJacobiKernel.py
      OUT_FILES CudaJacobiKernel2D.cu CudaJacobiKernel2D.h
      CudaJacobiKernel3D.cu CudaJacobiKernel3D.h)
waLBerla_compile_test( FILES codegen/CodegenJacobiGPU.cpp
                       DEPENDS blockforest timeloop gui CodegenJacobiGPUGeneratedCudaJacobiKernel )
waLBerla_execute_test( NAME CodegenJacobiGPU )

waLBerla_generate_target_from_python(NAME CodegenPoissonGPUGeneratedKernel FILE codegen/CudaPoisson.py
        OUT_FILES PoissonGPU.cu PoissonGPU.h )
waLBerla_compile_test( FILES codegen/CodegenPoissonGPU.cpp DEPENDS gui cuda timeloop CodegenPoissonGPUGeneratedKernel)
waLBerla_execute_test( NAME CodegenPoissonGPU )

# The following tests work only for CUDA enabled MPI
waLBerla_compile_test( FILES communication/CommTest )
#waLBerla_execute_test( NAME  CommTest PROCESSES 2)

waLBerla_compile_test( FILES CudaMPI DEPENDS blockforest timeloop gui )
#waLBerla_execute_test( NAME  CudaMPI )

waLBerla_compile_test( FILES AlignmentTest.cpp DEPENDS blockforest timeloop )

waLBerla_generate_target_from_python(NAME MicroBenchmarkGpuLbmGenerated FILE codegen/MicroBenchmarkGpuLbm.py
      OUT_FILES MicroBenchmarkStreamKernel.cu MicroBenchmarkCopyKernel.cu MicroBenchmarkStreamKernel.h MicroBenchmarkCopyKernel.h)
waLBerla_compile_test( FILES codegen/MicroBenchmarkGpuLbm.cpp DEPENDS MicroBenchmarkGpuLbmGenerated)

endif()