from pystencils_walberla.sweep import Sweep
from pystencils_walberla.cmake_integration import codegen

def jacobi2D(sweep):
    src = sweep.field("f1")
    dst = sweep.temporaryField(src)

    dst[0, 0] @= (src[1, 0] + src[-1, 0] + src[0, 1] + src[0, -1]) / (4 * S.h ** 2)

def jacobi3D(sweep):
    src = sweep.field("f1")
    dst = sweep.temporaryField(src)

    dst[0,0,0] @= (src[1,0,0] + src[-1,0,0] + src[0,1,0] + src[0, -1, 0] + src[0, 0, 1] + src[0, 0 , -1] ) / (6 * S.h**2)

Sweep.generate('CudaJacobiKernel2D', jacobi2D, dim=2, target='gpu')
Sweep.generate('CudaJacobiKernel3D', jacobi3D, dim=3, target='gpu')