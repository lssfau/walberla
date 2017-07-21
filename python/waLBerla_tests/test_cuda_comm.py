from waLBerla import *
import numpy as np
import pycuda.autoinit
from pycuda.gpuarray import *
from pycuda import *
from pystencils.field import createNumpyArrayWithLayout, getLayoutOfArray

blocks = createUniformBlockGrid( cells=(1,1,1), periodic=(1,1,1) )
cuda.addGpuFieldToStorage(blocks, "gpuField", float, fSize=1, ghostLayers=1, layout=field.fzyx, usePitchedMem=False)

gpuArr = cuda.toGpuArray(blocks[0]['gpuField'])

testField = createNumpyArrayWithLayout(gpuArr.shape, getLayoutOfArray(gpuArr))
testField[...] = 0
testField[1,1,1,0] = 1
gpuArr.set(testField)

scheme = createUniformBufferedScheme(blocks, "D3Q27")
scheme.addDataToCommunicate( cuda.createPackInfo(blocks, "gpuField") )

scheme()

gpuArr = cuda.toGpuArray(blocks[0]['gpuField'])

assert(np.allclose(np.ones([3,3,3,1]), gpuArr.get()))
