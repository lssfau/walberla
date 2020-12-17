import waLBerla
from waLBerla.field import toArray
import numpy as np


@waLBerla.callback("theCallback")
def theCallback(blocks):
    for block in blocks:
        np.copyto(toArray(block['srcIntFieldID']), toArray(block['dstIntFieldID']))
        np.copyto(toArray(block['srcDoubleFieldID']), toArray(block['dstDoubleFieldID']))
