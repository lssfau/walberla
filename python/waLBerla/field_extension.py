import numpy

# try:
#     from . import walberla_cpp
# except ImportError:
#     import walberla_cpp


# ----------------------------- Python functions to extend the C++ field module ---------------------------------

def normalizeGhostlayerInfo(field, withGhostLayers):
    """Takes one ghost layer parameter and returns an integer:
        True -> all ghost layers, False->no ghost layers"""

    def normalizeComponent(gl):
        if gl is False:
            return 0
        if gl is True:
            return field.nrOfGhostLayers
        if gl > field.nrOfGhostLayers:
            raise ValueError("Field only has %d ghost layers (requested %d)" % (field.nrOfGhostLayers, gl))
        return gl

    if hasattr(withGhostLayers, "__len__") and len(withGhostLayers) == 3:
        ghostLayers = [normalizeComponent(gl) for gl in withGhostLayers]
    else:
        ghostLayers = [normalizeComponent(withGhostLayers)] * 3
    return ghostLayers


def npArrayFromWaLBerlaField(field, withGhostLayers=False):
    """ Creates a numpy array view on the waLBerla field data
        @field: the waLBerla field
        @withGhostLayers: Possible values:
                            1. Boolean: False: no ghost layers included
                                        True:  all ghost layers included
                            2. Integer: number of ghost layers to include
                            3. List with three booleans or integers with ghost layer info for x,y,z direction
    """

    if not field:
        return None

    ghostLayers = normalizeGhostlayerInfo(field, withGhostLayers)

    if not hasattr(field, 'buffer'):  # Field adaptor -> create field with adapted values
        field = field.copyToField()

    if ghostLayers[0] == 0 and ghostLayers[1] == 0 and ghostLayers[2] == 0:
        return numpy.asarray(field.buffer(False))
    else:
        result = numpy.asarray(field.buffer(True))
        glCutoff = [field.nrOfGhostLayers - gl for gl in ghostLayers]
        view = result[glCutoff[0]:-glCutoff[0] if glCutoff[0] > 0 else None,
                      glCutoff[1]:-glCutoff[1] if glCutoff[1] > 0 else None,
                      glCutoff[2]:-glCutoff[2] if glCutoff[2] > 0 else None,
                      :]
        return view


def arrayFromWaLBerlaAdaptor(field, withGhostLayers=False):
    return npArrayFromWaLBerlaField(field.copyToField(), withGhostLayers)


def copyArrayToField(dstField, srcArray, slice=[slice(None, None, None)] * 3, withGhostLayers=False):
    """ Copies a numpy array into (part of) a waLBerla field

    Usually no copying has to take place between waLBerla fields and numpy arrays, since an array view can be
    constructed on a field that uses the same memory.
    When running certain numpy operations that cannot be done in-place,however, the data has to be copied back.

    @param dstField: waLBerla field, where the data is copied to
    @param srcArray: numpy array where to copy from
    @param slice:    the numpy array is allowed to be smaller than the field. In this case the target region
                     has to be specified via this 3 dimensional slice
    @param withGhostLayers: if true the ghost layers of the field are considered as well
    """
    dstAsArray = npArrayFromWaLBerlaField(dstField, withGhostLayers)
    numpy.copyto(dstAsArray[slice], srcArray)


def extend(cppFieldModule):
    def gatherField(blocks, blockDataName, sliceObj, allGather=False):
        field = cppFieldModule.gather(blocks, blockDataName, sliceObj, targetRank=-1 if allGather else 0)
        if field is not None:
            field = npArrayFromWaLBerlaField(field)
            field.flags.writeable = False
            return field
        else:
            return None

    cppFieldModule.toArray = npArrayFromWaLBerlaField
    cppFieldModule.adaptorToArray = arrayFromWaLBerlaAdaptor
    cppFieldModule.copyArrayToField = copyArrayToField
    cppFieldModule.gatherField = gatherField
