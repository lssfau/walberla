import sympy as sp
from lbmpy_walberla import Field, generateLatticeModelFiles, RefinementScaling

forceField = Field.createGeneric('force', spatialDimensions=3, indexDimensions=1, layout='c')
force = [forceField(0), forceField(1), forceField(2)]

omega = sp.Symbol("omega")

scaling = RefinementScaling()
scaling.addStandardRelaxationRateScaling(omega)
scaling.addForceScaling(forceField)

generateLatticeModelFiles(method='srt', stencil='D3Q19', forceModel='guo', force=force,
                          relaxationRates=[omega], refinementScaling=scaling)

