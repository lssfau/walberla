import waLBerla
from waLBerla.field import toArray
import numpy as np


@waLBerla.callback("theCallback")
def theCallback(blocks):
    for block in blocks:
        np.copyto(toArray(block['vec2Field']), toArray(block['sca2Field']))
        np.copyto(toArray(block['vec3Field']), toArray(block['sca3Field']))
