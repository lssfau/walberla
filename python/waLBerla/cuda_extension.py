from pycuda.gpuarray import GPUArray
import numpy as np
from .field_extension import normalizeGhostlayerInfo


def toGpuArray(f, withGhostLayers=True):
    """Converts a waLBerla GPUField to a pycuda GPUArray"""
    if not f:
        return None
    dtype = np.dtype(f.dtypeStr)
    strides = [dtype.itemsize * a for a in f.strides]
    res = GPUArray(f.sizeWithGhostLayers, dtype, gpudata=f.ptr, strides=strides)
    if withGhostLayers is True:
        return res

    ghostLayers = normalizeGhostlayerInfo(f, withGhostLayers)
    glCutoff = [f.nrOfGhostLayers - gl for gl in ghostLayers]
    res = res[glCutoff[0]:-glCutoff[0] if glCutoff[0] > 0 else None,
              glCutoff[1]:-glCutoff[1] if glCutoff[1] > 0 else None,
              glCutoff[2]:-glCutoff[2] if glCutoff[2] > 0 else None,
              :]
    return res


def extend(cppCudaModule):
    cppCudaModule.toGpuArray = toGpuArray
