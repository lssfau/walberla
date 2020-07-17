import numpy as np
import scipy
import scipy.ndimage

try:
    from .walberla_cpp import CellInterval, field
except ImportError:
    from walberla_cpp import CellInterval, field

from .core_extension import normalizeSlice


def setBoundaryFromArray(blocks, boundaryID, targetSlice, imageArr, boundaryConfig,
                         resizeFunc=None, extrusionCoordinate=-1):
    """Initializes Boundary Handling using an image
    :param blocks: the block storage
    :param boundaryID: block data name of boundary handling
    :param targetSlice: slice of the domain where the image should be placed. If a 2D slice is given, the image is
                        automatically extruded along the third coordinate. For 3D (or 1D) slices an extrusionCoordinate
                        has to be given.
                        Example: for targetSlice=[0.25:0.75,  0  , 0.25:0.75] the image is resized to half
                        the x-z domain size, placed in the middle of the domain and extruded in y direction.
    :param imageArr:    a 2D array used to set up the boundaries
    :param boundaryConfig: dictionary mapping values of imageArr to boundary configurations,
                           used as index array in forceBoundary()
    :param resizeFunc: if the given slice does not match the shape of imageArr, the image array has to be resized.
                       This can not be done automatically since interpolation would change the values, and the mapping
                       given in boundaryConfig is incorrect. Thus the resize function has to be supplied by the caller.
                       Here a function "resize(imgArr, newSize)" has to be passed, that returns a resized array of
                       shape "newSize".
    :param extrusionCoordinate: only necessary if 3D or 1D slices are given for targetSlice.
                                See documentation of targetSlice.
    """

    if len(blocks) == 0:
        return

    nrOfGhostLayers = blocks[0][boundaryID].getFlagField().nrOfGhostLayers

    imageArr = np.rot90(imageArr, 3)

    sliceWithGhostLayers = 'g' in targetSlice
    targetSlice = [s for s in targetSlice if s != 'g']
    size = [s + 2 * nrOfGhostLayers if sliceWithGhostLayers else s for s in blocks.getDomainCellBB().size]
    imageCellInterval = CellInterval.fromSlice(normalizeSlice(targetSlice, size))
    if sliceWithGhostLayers:
        imageCellInterval.shift(-nrOfGhostLayers, -nrOfGhostLayers, -nrOfGhostLayers)

    # Automatic detection of extrusion coordinate
    if extrusionCoordinate < 0 or extrusionCoordinate > 2:
        possibleExtrusionCoordinate = np.array([0, 0, 0])
        for i in range(3):
            if imageCellInterval.min[i] == imageCellInterval.max[i]:
                possibleExtrusionCoordinate[i] = 1
                extrusionCoordinate = i
        if sum(possibleExtrusionCoordinate) != 1:
            raise ValueError("No valid extrusionCoordinate given - "
                             "and extrusion coordinate could not be found automatically")
    assert (extrusionCoordinate < 3 and extrusionCoordinate >= 0)

    # Resize image
    imageBounds = list(imageCellInterval.size)
    del imageBounds[extrusionCoordinate]
    if imageArr.shape != tuple(imageBounds):
        if resizeFunc is None:
            raise ValueError("The given image size does not match the target slice: "
                             "resizing would be necessary but no resizeFunc was given.")
        imageArr = resizeFunc(imageArr, imageBounds)

    assert imageArr.shape == tuple(imageBounds)
    assert imageArr.dtype.kind == 'i', "imageArr has to be of integer type"

    unusedIdx = 0
    while unusedIdx in boundaryConfig:
        unusedIdx += 1

    def make2Dfrom3DSlice(targetSlice, extrusionCoordinate):
        l = list(targetSlice)
        del l[extrusionCoordinate]
        return l

    for block in blocks:
        blockCellInterval = blocks.getBlockCellBB(block)
        blockCellInterval.expand(nrOfGhostLayers)
        intersectionGlobalCoord = blockCellInterval.getIntersection(imageCellInterval)

        if intersectionGlobalCoord.empty():
            continue
        intersectionLocalCoord = blocks.transformGlobalToLocal(block, intersectionGlobalCoord)

        # Create a field with same size as block
        blockCellBB = blocks.getBlockCellBB(block)
        wlbIndexField = field.createField(list(blockCellBB.size), np.int32, ghostLayers=nrOfGhostLayers)
        indexField = field.toArray(wlbIndexField, withGhostLayers=nrOfGhostLayers)[:, :, :, :]
        indexField[:, :, :, :] = unusedIdx

        # Copy image into this field
        targetSlice = intersectionLocalCoord.getShifted(nrOfGhostLayers, nrOfGhostLayers, nrOfGhostLayers).toSlice()

        minCoord = np.array(imageCellInterval.min)
        imgTargetSlice = intersectionGlobalCoord.getShifted(*(-minCoord)).toSlice()
        sliceInImage = make2Dfrom3DSlice(imgTargetSlice, extrusionCoordinate)
        indexField[targetSlice + [0]] = imageArr[sliceInImage]

        block[boundaryID].forceBoundary(wlbIndexField, boundaryConfig)


def binaryResize(img, newSize):
    """This can be used as resize function for setBoundaryFromArray for arrays with
       zero and ones. After resizing the image is again binarized"""
    img = scipy.misc.imresize(img, size=newSize)
    img[img <= 254] = 0
    img[img > 254] = 1
    img = img.astype(np.int32)
    return img


def setBoundaryFromBlackAndWhiteImage(blocks, boundaryID, targetSlice, imagePath, boundaryConfig,
                                      extrusionCoordinate=-1):
    """Loads array from image file and calls setBoundaryFromArray.

     :param imagePath: path to image file.

     For the other parameters see documentation of setBoundaryFromArray.
    """
    imgArr = scipy.ndimage.imread(imagePath, flatten=True).astype(int)
    setBoundaryFromArray(blocks, boundaryID, targetSlice, imgArr, {0: boundaryConfig},
                         resizeFunc=binaryResize, extrusionCoordinate=extrusionCoordinate)


def setFieldUsingFlagMask(blocks, targetField, targetValue, flagField, flagNames):
    """
    Sets all values of a target field to given value where a certain flag is set.

    :param blocks:  the block structure
    :param targetField: the field that is modified
    :param targetValue: value that is written to all entries of cells where the given flag is set
    :param flagField: the flag field
    :param flagNames: list of flag names. If one of the flags is set, the target value is written
    """
    for b in blocks:
        mask = 0
        for flagName in flagNames:
            mask |= b[flagField].flag(flagName)

        targetArr = field.toArray(b[targetField], True)
        flagArr = field.toArray(b[flagField], True)[:, :, :, 0]
        targetArr[np.bitwise_and(flagArr, mask) > 0, :] = targetValue
