###################################################################################################
#
# Module gpu
#
###################################################################################################

add_library( gpu )
target_link_libraries( gpu PUBLIC blockforest core communication domain_decomposition executiontree field stencil lbm )

# sources for HIP and CUDA
target_sources( gpu
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
      )

# sources only for CUDA
if (WALBERLA_BUILD_WITH_CUDA)
target_sources( gpu
        PRIVATE
        NVTX.h
        )
endif (WALBERLA_BUILD_WITH_CUDA)

add_subdirectory( sweeps )
add_subdirectory( communication )
add_subdirectory( lbm )

###################################################################################################