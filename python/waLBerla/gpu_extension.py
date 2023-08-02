import cupy as cp
from cupy.cuda import MemoryPointer, UnownedMemory
import numpy as np
from .field_extension import normalize_ghostlayer_info


def to_gpu_array(f, with_ghost_layers=True):
    """Converts a waLBerla GPUField to a cupy ndarray"""
    if not f:
        return None
    dtype = np.dtype(f.dtypeStr)
    strides = [dtype.itemsize * a for a in f.strides]

    allocated_bytes = np.prod(f.allocSize) * dtype.itemsize
    memory_pointer = MemoryPointer(UnownedMemory(f.ptr, allocated_bytes, f), 0)

    res = cp.ndarray(shape=f.sizeWithGhostLayers, dtype=dtype,
                     memptr=memory_pointer, strides=strides)

    if with_ghost_layers is True:
        return res

    ghost_layers = normalize_ghostlayer_info(f, with_ghost_layers)
    cutoff = [f.nrOfGhostLayers - gl for gl in ghost_layers]
    res = res[slice(cutoff[0], -cutoff[0], 1) if cutoff[0] > 0 else slice(None, None, None),
              slice(cutoff[1], -cutoff[1], 1) if cutoff[1] > 0 else slice(None, None, None),
              slice(cutoff[2], -cutoff[2], 1) if cutoff[2] > 0 else slice(None, None, None),
              slice(None, None, None)]
    return res


def extend(cpp_gpu_module):
    cpp_gpu_module.toGpuArray = to_gpu_array
