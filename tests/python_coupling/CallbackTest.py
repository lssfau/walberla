import waLBerla
# import numpy as np


@waLBerla.callback("cb1")
def someCallback(input1, input2):
    return input1 + input2


@waLBerla.callback("cb2")
def fieldCallback(field):
    npArray = waLBerla.field.toArray(field)
    npArray[0, 0, 0] = 42

    npArrayGl = waLBerla.field.toArray(field, withGhostLayers=True)
    print(npArrayGl.shape)
    npArrayGl[0, 0, 0] = 5
