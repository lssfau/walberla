###################################################################################################
#
# Tests for gpu
#
###################################################################################################

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_HIP )

   waLBerla_compile_test( FILES communication/GPUPackInfoTest.cpp DEPENDS walberla::blockforest )
waLBerla_execute_test( NAME  GPUPackInfoTest )

   waLBerla_compile_test( FILES communication/GPUPackInfoCommunicationTest.cpp
         DEPENDS walberla::domain_decomposition walberla::blockforest walberla::stencil )
waLBerla_execute_test( NAME  GPUPackInfoCommunicationTest )

   waLBerla_compile_test( FILES communication/GPUBlockSelectorCommunicationTest.cpp
         DEPENDS walberla::domain_decomposition walberla::blockforest walberla::stencil )
waLBerla_execute_test( NAME  GPUBlockSelectorCommunicationTest )

waLBerla_compile_test(
   NAME TestGenericGpuPackInfos
   FILES communication/TestGenericGpuPackInfos.cpp
   DEPENDS walberla::blockforest walberla::stencil walberla::gpu
)

waLBerla_execute_test( NAME TestGenericGpuPackInfos_1 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 1 )
waLBerla_execute_test( NAME TestGenericGpuPackInfos_2 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 2 )
waLBerla_execute_test( NAME TestGenericGpuPackInfos_4 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 4 )


waLBerla_compile_test( FILES FieldTransferTest.cpp )
waLBerla_execute_test( NAME  FieldTransferTest )

if( WALBERLA_BUILD_WITH_CUDA )
   waLBerla_compile_test( FILES TestGpuRAII.cu DEPENDS walberla::gpu )
elseif( WALBERLA_BUILD_WITH_HIP )
waLBerla_compile_test( FILES TestGpuRAII.hip DEPENDS walberla::gpu )
endif()

waLBerla_execute_test( NAME TestGpuRAII )

if( WALBERLA_BUILD_WITH_CUDA)
# These tests include handwritten CUDA code and won't compile with HIP

   waLBerla_compile_test( FILES SimpleKernelTest.cpp Kernels.cu DEPENDS walberla::blockforest walberla::timeloop )
waLBerla_execute_test( NAME  SimpleKernelTest )

waLBerla_compile_test( FILES FieldIndexing3DTest.cpp FieldIndexing3DTest.cu )
waLBerla_execute_test( NAME  FieldIndexing3DTest )


if( WALBERLA_BUILD_WITH_CODEGEN )
waLBerla_generate_target_from_python(NAME CodegenJacobiGPUGeneratedCudaJacobiKernel FILE codegen/CudaJacobiKernel.py
      OUT_FILES CudaJacobiKernel2D.cu CudaJacobiKernel2D.h
      CudaJacobiKernel3D.cu CudaJacobiKernel3D.h)
waLBerla_compile_test( FILES codegen/CodegenJacobiGPU.cpp
      DEPENDS walberla::blockforest walberla::timeloop CodegenJacobiGPUGeneratedCudaJacobiKernel )
waLBerla_execute_test( NAME CodegenJacobiGPU )

waLBerla_generate_target_from_python(NAME CodegenPoissonGPUGeneratedKernel FILE codegen/CudaPoisson.py
        OUT_FILES PoissonGPU.cu PoissonGPU.h )
   waLBerla_compile_test( FILES codegen/CodegenPoissonGPU.cpp DEPENDS walberla::gpu walberla::timeloop CodegenPoissonGPUGeneratedKernel )
waLBerla_execute_test( NAME CodegenPoissonGPU )
endif()

# The following tests work only for CUDA enabled MPI
waLBerla_compile_test( FILES communication/CommTest.cpp )
#waLBerla_execute_test( NAME  CommTest PROCESSES 2)

   waLBerla_compile_test( FILES CudaMPI.cpp DEPENDS walberla::blockforest walberla::timeloop )
#waLBerla_execute_test( NAME  CudaMPI )

endif()

   waLBerla_compile_test( FILES AlignmentTest.cpp DEPENDS walberla::blockforest walberla::timeloop )

if( WALBERLA_BUILD_WITH_CODEGEN )
waLBerla_generate_target_from_python(NAME MicroBenchmarkGpuLbmGenerated FILE codegen/MicroBenchmarkGpuLbm.py
      OUT_FILES MicroBenchmarkStreamKernel.${CODEGEN_FILE_SUFFIX} MicroBenchmarkCopyKernel.${CODEGEN_FILE_SUFFIX} MicroBenchmarkStreamKernel.h MicroBenchmarkCopyKernel.h)
waLBerla_compile_test( FILES codegen/MicroBenchmarkGpuLbm.cpp DEPENDS MicroBenchmarkGpuLbmGenerated)

waLBerla_generate_target_from_python(NAME CodegenGeneratedGPUFieldPackInfo FILE codegen/GeneratedFieldPackInfoTestGPU.py
        OUT_FILES ScalarFieldCommunicationGPU.${CODEGEN_FILE_SUFFIX} ScalarFieldCommunicationGPU.h
        ScalarFieldPullReductionGPU.${CODEGEN_FILE_SUFFIX} ScalarFieldPullReductionGPU.h )
waLBerla_compile_test( FILES codegen/GeneratedFieldPackInfoTestGPU.cpp
      DEPENDS walberla::blockforest walberla::core walberla::field CodegenGeneratedGPUFieldPackInfo )
waLBerla_execute_test( NAME GeneratedFieldPackInfoTestGPU )


waLBerla_generate_target_from_python(NAME CodegenGeneratedGPUFieldNonUniformPackInfo FILE codegen/GeneratedFieldNonUniformPackInfoTestGPU.py
        OUT_FILES ScalarFieldNonUniformCommunicationGPU.${CODEGEN_FILE_SUFFIX} ScalarFieldNonUniformCommunicationGPU.h)
waLBerla_compile_test( FILES codegen/GeneratedFieldNonUniformPackInfoTestGPU.cpp
        DEPENDS walberla::blockforest walberla::core walberla::field CodegenGeneratedGPUFieldNonUniformPackInfo )
waLBerla_execute_test( NAME GeneratedFieldNonUniformPackInfoTestGPU )
endif()


    if (WALBERLA_BUILD_WITH_PYTHON)

        waLBerla_link_files_to_builddir( *.py )

        waLBerla_compile_test( FILES TestShiftedPeriodicityGPU.cpp DEPENDS walberla::gpu walberla::blockforest walberla::field walberla::python_coupling )
        waLBerla_execute_test( NAME TestShiftedPeriodicityGPU COMMAND $<TARGET_FILE:TestShiftedPeriodicityGPU> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetupGPU.py )

    endif()

endif()
