import sympy as sp
from lbmpy.boundaries import NoSlip, UBB
from lbmpy_walberla import Field, generateLatticeModelFiles, RefinementScaling
from lbmpy.creationfunctions import createLatticeBoltzmannMethod
from lbmpy_walberla.boundary import createBoundaryClass
from pystencils_walberla.cmake_integration import codegen

# ------------- Lattice Model ------------------------------
forceField = Field.createGeneric('force', spatialDimensions=3, indexDimensions=1, layout='fzyx')
force = [forceField(0), forceField(1), forceField(2)]

omega = sp.Symbol("omega")

scaling = RefinementScaling()
scaling.addStandardRelaxationRateScaling(omega)
scaling.addForceScaling(forceField)

generateLatticeModelFiles(className='SrtWithForceFieldModel',
                          method='srt', stencil='D3Q19', forceModel='guo', force=force,
                          relaxationRates=[omega], refinementScaling=scaling)


def genBoundary():
    boundary = UBB([0.05, 0, 0], dim=3, name="MyUBB")
    method = createLatticeBoltzmannMethod(stencil='D3Q19', method='srt')
    return createBoundaryClass(boundary, method)

def genNoSlip():
    boundary = NoSlip(name='MyNoSlip')
    method = createLatticeBoltzmannMethod(stencil='D3Q19', method='srt')
    return createBoundaryClass(boundary, method)

codegen.register(['MyUBB.h', 'MyUBB.cpp'], genBoundary)
codegen.register(['MyNoSlip.h', 'MyNoSlip.cpp',], genNoSlip)

