from pycuda.gpuarray import GPUArray
import numpy as np
from .field_extension import normalize_ghostlayer_info


def to_gpu_array(f, with_ghost_layers=True):
    """Converts a waLBerla GPUField to a pycuda GPUArray"""
    if not f:
        return None
    dtype = np.dtype(f.dtypeStr)
    strides = [dtype.itemsize * a for a in f.strides]
    res = GPUArray(f.sizeWithGhostLayers, dtype, gpudata=f.ptr, strides=strides)
    if with_ghost_layers is True:
        return res

    ghost_layers = normalize_ghostlayer_info(f, with_ghost_layers)
    cutoff = [f.nrOfGhostLayers - gl for gl in ghost_layers]
    res = res[cutoff[0]:-cutoff[0] if cutoff[0] > 0 else None,
              cutoff[1]:-cutoff[1] if cutoff[1] > 0 else None,
              cutoff[2]:-cutoff[2] if cutoff[2] > 0 else None,
              :]
    return res


def extend(cpp_cuda_module):
    cpp_cuda_module.toGpuArray = to_gpu_array
