###################################################################################################
#
# Module gpu
#
###################################################################################################

add_library( walberla_gpu )
add_library( walberla::gpu ALIAS walberla_gpu )
target_link_libraries( walberla_gpu PUBLIC walberla::blockforest walberla::core walberla::communication walberla::domain_decomposition walberla::executiontree walberla::field walberla::stencil walberla::lbm )

if ( WALBERLA_BUILD_WITH_CUDA )
   target_link_libraries( walberla_gpu PUBLIC CUDA::cudart)
endif ( WALBERLA_BUILD_WITH_CUDA )

if ( WALBERLA_BUILD_WITH_HIP )
    target_link_libraries( walberla_gpu PUBLIC hip::host )
endif ( WALBERLA_BUILD_WITH_HIP )

# sources for HIP and CUDA
target_sources( walberla_gpu
      PRIVATE
      AlignedAllocation.h
      AddGPUFieldToStorage.h
      ErrorChecking.h
      FieldCopy.h
      GPUCopy.cpp
      FieldIndexingXYZ.h
      FieldIndexing3D.h
      AddGPUFieldToStorage.impl.h
      GPUField.h
      GPUWrapper.h
      DeviceWrapper.h
      FieldAccessor3D.h
      DeviceSelectMPI.h
      HostFieldAllocator.h
      FieldAccessor.h
      FieldIndexing3D.impl.h
      GPUCopy.h
      FieldAccessorXYZ.h
      FieldIndexingXYZ.impl.h
      FieldIndexing.h
      AlignedAllocation.cpp
      GPUField.impl.h
      ParallelStreams.cpp
      FieldIndexing.impl.h
      Kernel.h
      ParallelStreams.h
      GPURAII.h
      DeviceSelectMPI.cpp
      ShiftedPeriodicity.cu
      ShiftedPeriodicity.h
      )

# sources only for CUDA
if (WALBERLA_BUILD_WITH_CUDA)
   target_sources( walberla_gpu
        PRIVATE
        NVTX.h
        )
endif (WALBERLA_BUILD_WITH_CUDA)

add_subdirectory( sweeps )
add_subdirectory( communication )
add_subdirectory( lbm )

###################################################################################################