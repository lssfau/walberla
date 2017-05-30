from pycuda.gpuarray import GPUArray
import numpy as np

def toGpuArray(f):
    """Converts a waLBerla GPUField to a pycuda GPUArray"""
    if not f:
        return None
    dtype = np.dtype(f.dtypeStr)
    strides = [dtype.itemsize*a for a in f.strides]
    return GPUArray(f.sizeWithGhostLayers, dtype, gpudata=f.ptr, strides=strides)


def extend(cppCudaModule):
    cppCudaModule.toGpuArray = toGpuArray

