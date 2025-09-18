###################################################################################################
#
# Tests for gpu
#
###################################################################################################

if( WALBERLA_BUILD_WITH_CUDA OR WALBERLA_BUILD_WITH_HIP )

    waLBerla_add_test_executable ( GPUPackInfoTest communication/GPUPackInfoTest.cpp)
    target_link_libraries ( GPUPackInfoTest walberla::blockforest walberla::core walberla::gpu )
    waLBerla_execute_test( NAME  GPUPackInfoTest )

    waLBerla_add_test_executable ( GPUPackInfoCommunicationTest communication/GPUPackInfoCommunicationTest.cpp)
    target_link_libraries ( GPUPackInfoCommunicationTest walberla::domain_decomposition walberla::blockforest
            walberla::stencil walberla::core walberla::gpu )
    waLBerla_execute_test( NAME  GPUPackInfoCommunicationTest )

    waLBerla_add_test_executable ( GPUBlockSelectorCommunicationTest communication/GPUBlockSelectorCommunicationTest.cpp)
    target_link_libraries ( GPUBlockSelectorCommunicationTest walberla::domain_decomposition walberla::blockforest
            walberla::stencil walberla::core walberla::gpu )
    waLBerla_execute_test( NAME  GPUBlockSelectorCommunicationTest )

    waLBerla_add_test_executable ( TestGenericGpuPackInfos communication/TestGenericGpuPackInfos.cpp )
    target_link_libraries ( TestGenericGpuPackInfos walberla::blockforest walberla::stencil walberla::gpu walberla::core )

    waLBerla_execute_test( NAME TestGenericGpuPackInfos_1 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 1 )
    waLBerla_execute_test( NAME TestGenericGpuPackInfos_2 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 2 )
    waLBerla_execute_test( NAME TestGenericGpuPackInfos_4 COMMAND $<TARGET_FILE:TestGenericGpuPackInfos> PROCESSES 4 )


    waLBerla_add_test_executable ( FieldTransferTest FieldTransferTest.cpp )
    target_link_libraries ( FieldTransferTest walberla::gpu walberla::core )
    waLBerla_execute_test( NAME  FieldTransferTest )

    if( WALBERLA_BUILD_WITH_CUDA )
        waLBerla_add_test_executable ( TestGpuRAII TestGpuRAII.cu )
        target_link_libraries ( TestGpuRAII walberla::gpu walberla::core )
    elseif( WALBERLA_BUILD_WITH_HIP )
        waLBerla_add_test_executable ( TestGpuRAII TestGpuRAII.hip )
        target_link_libraries ( TestGpuRAII walberla::gpu walberla::core )
    endif()

    waLBerla_execute_test( NAME TestGpuRAII )

    if( WALBERLA_BUILD_WITH_CUDA)
        # These tests include handwritten CUDA code and won't compile with HIP

        waLBerla_add_test_executable ( SimpleKernelTest SimpleKernelTest.cpp Kernels.cu )
        target_link_libraries( SimpleKernelTest walberla::blockforest walberla::timeloop walberla::core walberla::gpu )
        waLBerla_execute_test( NAME  SimpleKernelTest )

        waLBerla_add_test_executable ( FieldIndexing3DTest FieldIndexing3DTest.cpp FieldIndexing3DTest.cu )
        target_link_libraries( FieldIndexing3DTest walberla::core walberla::gpu )
        waLBerla_execute_test( NAME  FieldIndexing3DTest )


        if ( WALBERLA_BUILD_WITH_CODEGEN )
            waLBerla_generate_target_from_python(NAME CodegenJacobiGPUGeneratedCudaJacobiKernel
                    FILE codegen/CudaJacobiKernel.py
                    OUT_FILES CudaJacobiKernel2D.cu CudaJacobiKernel2D.h CudaJacobiKernel3D.cu CudaJacobiKernel3D.h
            )

            waLBerla_add_test_executable ( CodegenJacobiGPU codegen/CodegenJacobiGPU.cpp )
            target_link_libraries ( CodegenJacobiGPU walberla::blockforest walberla::timeloop
                    CodegenJacobiGPUGeneratedCudaJacobiKernel walberla::core walberla::gpu )
            waLBerla_execute_test( NAME CodegenJacobiGPU )

            waLBerla_generate_target_from_python(NAME CodegenPoissonGPUGeneratedKernel FILE codegen/CudaPoisson.py
                OUT_FILES PoissonGPU.cu PoissonGPU.h
            )

            waLBerla_add_test_executable ( CodegenPoissonGPU codegen/CodegenPoissonGPU.cpp )
            target_link_libraries ( CodegenPoissonGPU walberla::gpu walberla::timeloop CodegenPoissonGPUGeneratedKernel
                    walberla::core walberla::gpu )
            waLBerla_execute_test( NAME CodegenPoissonGPU )
        endif()

        # The following tests work only for CUDA enabled MPI
        waLBerla_add_test_executable ( CommTest communication/CommTest.cpp )
        target_link_libraries ( CommTest walberla::core walberla::gpu )
        #waLBerla_execute_test( NAME  CommTest PROCESSES 2)

        waLBerla_add_test_executable ( CudaMPI CudaMPI.cpp )
        target_link_libraries ( CudaMPI walberla::blockforest walberla::timeloop walberla::core walberla::gpu )
        #waLBerla_execute_test( NAME  CudaMPI )

    endif()

    waLBerla_add_test_executable ( AlignmentTest AlignmentTest.cpp )
    target_link_libraries ( AlignmentTest walberla::blockforest walberla::timeloop walberla::core walberla::gpu )

    if( WALBERLA_BUILD_WITH_CODEGEN )
        waLBerla_generate_target_from_python(NAME MicroBenchmarkGpuLbmGenerated
                FILE codegen/MicroBenchmarkGpuLbm.py
                OUT_FILES MicroBenchmarkStreamKernel.${CODEGEN_FILE_SUFFIX}
                MicroBenchmarkCopyKernel.${CODEGEN_FILE_SUFFIX} MicroBenchmarkStreamKernel.h MicroBenchmarkCopyKernel.h
        )

        waLBerla_add_test_executable ( MicroBenchmarkGpuLbm codegen/MicroBenchmarkGpuLbm.cpp )
        target_link_libraries ( MicroBenchmarkGpuLbm MicroBenchmarkGpuLbmGenerated walberla::core walberla::gpu )

        waLBerla_generate_target_from_python(NAME CodegenGeneratedGPUFieldPackInfo
                FILE codegen/GeneratedFieldPackInfoTestGPU.py
                OUT_FILES ScalarFieldCommunicationGPU.${CODEGEN_FILE_SUFFIX} ScalarFieldCommunicationGPU.h
                ScalarFieldPullReductionGPU.${CODEGEN_FILE_SUFFIX} ScalarFieldPullReductionGPU.h
        )

        waLBerla_add_test_executable ( GeneratedFieldPackInfoTestGPU codegen/GeneratedFieldPackInfoTestGPU.cpp)
        target_link_libraries ( GeneratedFieldPackInfoTestGPU walberla::blockforest walberla::core walberla::field
                CodegenGeneratedGPUFieldPackInfo walberla::core walberla::gpu )
        waLBerla_execute_test( NAME GeneratedFieldPackInfoTestGPU )


        waLBerla_generate_target_from_python(NAME CodegenGeneratedGPUFieldNonUniformPackInfo
                FILE codegen/GeneratedFieldNonUniformPackInfoTestGPU.py
                OUT_FILES ScalarFieldNonUniformCommunicationGPU.${CODEGEN_FILE_SUFFIX}
                ScalarFieldNonUniformCommunicationGPU.h
        )

        waLBerla_add_test_executable ( GeneratedFieldNonUniformPackInfoTestGPU codegen/GeneratedFieldNonUniformPackInfoTestGPU.cpp)
        target_link_libraries ( GeneratedFieldNonUniformPackInfoTestGPU walberla::blockforest walberla::core
                walberla::field CodegenGeneratedGPUFieldNonUniformPackInfo walberla::core walberla::gpu )
        waLBerla_execute_test( NAME GeneratedFieldNonUniformPackInfoTestGPU )
endif()


    if (WALBERLA_BUILD_WITH_PYTHON)

        waLBerla_link_files_to_builddir( *.py )

        waLBerla_add_test_executable ( TestShiftedPeriodicityGPU TestShiftedPeriodicityGPU.cpp )
        target_link_libraries ( TestShiftedPeriodicityGPU walberla::gpu walberla::blockforest walberla::field
                walberla::python_coupling walberla::core walberla::gpu )
        waLBerla_execute_test( NAME TestShiftedPeriodicityGPU
                COMMAND $<TARGET_FILE:TestShiftedPeriodicityGPU> ${CMAKE_CURRENT_SOURCE_DIR}/TestShiftedPeriodicitySetupGPU.py )

    endif()

endif()
