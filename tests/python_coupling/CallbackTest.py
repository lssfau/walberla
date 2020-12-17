import waLBerla
# import numpy as np


@waLBerla.callback("cb1")
def someCallback(input1, input2):
    return input1 + input2


@waLBerla.callback("cb2")
def fieldCallback(field):
    numpy_array = waLBerla.field.toArray(field)
    numpy_array[0, 0, 0] = 42

    numpy_array_with_gl = waLBerla.field.toArray(field, with_ghost_layers=True)
    print(numpy_array_with_gl.shape)
    numpy_array_with_gl[0, 0, 0] = 5
