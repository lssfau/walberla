###################################################################################################
#
# Module gpu
#
###################################################################################################

add_library( walberla_gpu )
add_library( walberla::gpu ALIAS walberla_gpu )
target_link_libraries( walberla_gpu PUBLIC walberla::blockforest walberla::core walberla::communication
        walberla::domain_decomposition walberla::executiontree walberla::field walberla::stencil walberla::lbm )

# CUDA linkage
if ( WALBERLA_BUILD_WITH_CUDA )
    #include directories and cudart lib is needed for cpp files that use cuda headers/libs
    target_link_libraries ( walberla_gpu PUBLIC CUDA::cudart )
    target_link_libraries ( walberla_gpu PUBLIC CUDA::nvtx3 )
endif ( WALBERLA_BUILD_WITH_CUDA )

# HIP linkage
if ( WALBERLA_BUILD_WITH_HIP )
    target_link_libraries( walberla_gpu PUBLIC hip::host )
endif ( WALBERLA_BUILD_WITH_HIP )

# sources for HIP and CUDA
target_sources( walberla_gpu
        PUBLIC
        AlignedAllocation.h
        AddGPUFieldToStorage.h
        AddGPUFieldToStorage.impl.h
        ErrorChecking.h
        FieldCopy.h
        FieldIndexingXYZ.h
        FieldIndexingXYZ.impl.h
        FieldIndexing3D.h
        FieldIndexing3D.impl.h
        GPUField.h
        GPUField.impl.h
        GPUWrapper.h
        DeviceWrapper.h
        FieldAccessor3D.h
        DeviceSelectMPI.h
        HostFieldAllocator.h
        FieldAccessor.h
        GPUCopy.h
        FieldAccessorXYZ.h
        FieldIndexing.h
        FieldIndexing.impl.h
        Kernel.h
        ParallelStreams.h
        GPURAII.h
        ShiftedPeriodicity.h
)

target_sources( walberla_gpu
        PRIVATE
        GPUCopy.cpp
        AlignedAllocation.cpp
        ParallelStreams.cpp
        DeviceSelectMPI.cpp
        ShiftedPeriodicity.cu
)

# sources only for CUDA
if (WALBERLA_BUILD_WITH_CUDA)
    target_sources( walberla_gpu
            PUBLIC
            NVTX.h
    )
endif (WALBERLA_BUILD_WITH_CUDA)

add_subdirectory( sweeps )
add_subdirectory( communication )
add_subdirectory( lbm )

###################################################################################################