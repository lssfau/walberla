###################################################################################################
#
# Module cuda
#
###################################################################################################

add_library( cuda )
target_link_libraries( cuda PUBLIC blockforest core communication domain_decomposition executiontree field stencil lbm )
target_sources( cuda
      PRIVATE
      AlignedAllocation.h
      AddGPUFieldToStorage.h
      ErrorChecking.h
      ExecutionTreeGPU.h
      FieldCopy.h
      GPUCopy.cpp
      NVTX.h
      FieldIndexingXYZ.h
      FieldIndexing3D.h
      AddGPUFieldToStorage.impl.h
      GPUField.h
      FieldAccessor3D.h
      DeviceSelectMPI.h
      HostFieldAllocator.h
      FieldAccessor.h
      FieldIndexing3D.impl.h
      GPUCopy.h
      FieldAccessorXYZ.h
      FieldIndexingXYZ.impl.h
      ExecutionTreeSweepGPU.h
      FieldIndexing.h
      AlignedAllocation.cpp
      GPUField.impl.h
      ParallelStreams.cpp
      FieldIndexing.impl.h
      Kernel.h
      ParallelStreams.h
      CudaRAII.h
      DeviceSelectMPI.cpp
      )

add_subdirectory( sweeps )
add_subdirectory( communication )
add_subdirectory( lbm )

###################################################################################################