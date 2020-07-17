from waLBerla import field, createUniformBlockGrid, createUniformBufferedScheme, cuda
import numpy as np
import pycuda.autoinit  # noqa: F401
import pycuda.gpuarray as gpuArr
# from pycuda import *
from pystencils.field import createNumpyArrayWithLayout, getLayoutOfArray

blocks = createUniformBlockGrid(cells=(1, 1, 1), periodic=(1, 1, 1))
cuda.addGpuFieldToStorage(blocks, "gpuField", float, fSize=1, ghostLayers=1, layout=field.fzyx, usePitchedMem=False)

gpuArr = cuda.toGpuArray(blocks[0]['gpuField'])  # noqa: F811

testField = createNumpyArrayWithLayout(gpuArr.shape, getLayoutOfArray(gpuArr))
testField[...] = 0
testField[1, 1, 1, 0] = 1
gpuArr.set(testField)

scheme = createUniformBufferedScheme(blocks, "D3Q27")
scheme.addDataToCommunicate(cuda.createPackInfo(blocks, "gpuField"))

scheme()

gpuArr = cuda.toGpuArray(blocks[0]['gpuField'])

assert (np.allclose(np.ones([3, 3, 3, 1]), gpuArr.get()))
