from waLBerla import makeSlice, field, mpi, lbm, createUniformBlockGrid, createUniformBufferedScheme
from waLBerla.geometry_setup import setBoundaryFromBlackAndWhiteImage, setFieldUsingFlagMask
import itertools

import os
import numpy as np
import scipy

imageFile = os.path.join(os.path.dirname(__file__), 'wing.png')


def setBoundariesChannel(blocks, boundaryHandlingID):
    for block in blocks:
        b = block[boundaryHandlingID]
        if block.atDomainMinBorder[1]:
            b.forceBoundary('NoSlip', makeSlice[:, 0, :, 'g'])
        if block.atDomainMaxBorder[1]:
            b.forceBoundary('NoSlip', makeSlice[:, -1, :, 'g'])
        b.fillWithDomain()


class ForceCalculationMasks:
    @staticmethod
    def addToBlock(block, blockStorage):
        pdfFieldArr = field.toArray(block['pdfs'])
        flagFieldArr = field.toArray(block['flags'])[:, :, :, 0]
        directions = block['pdfs'].latticeModel.directions
        maskArr = np.zeros(pdfFieldArr.shape, dtype=bool)
        pdfDirectionArr = np.zeros(list(pdfFieldArr.shape) + [3])

        nearBoundaryFlag = block['flags'].flag("fluid")
        noSlipFlag = block['flags'].flag("NoSlip")

        innerPartOfDomain = itertools.product(range(2, maskArr.shape[0] - 2),
                                              range(2, maskArr.shape[1] - 2),
                                              range(maskArr.shape[2]))

        for x, y, z in innerPartOfDomain:
            if flagFieldArr[x, y, z] & nearBoundaryFlag:
                for dirIdx, dir in enumerate(directions):
                    nx, ny, nz = x + dir[0], y + dir[1], z + dir[2]
                    if flagFieldArr[nx, ny, nz] & noSlipFlag:
                        maskArr[x, y, z, dirIdx] = True
                        pdfDirectionArr[x, y, z, :] = dir
        return ForceCalculationMasks(maskArr, pdfDirectionArr)

    def __init__(self, maskArr, pdfDirectionArr):
        self._maskArr = maskArr
        self._pdfDirectionArr = pdfDirectionArr

    def calculateForceOnBoundary(self, pdfField):
        force = np.array([0.0] * 3)
        pdfFieldArr = field.toArray(pdfField)
        for i in range(3):
            fArr = pdfFieldArr[self._maskArr] * self._pdfDirectionArr[self._maskArr, i]
            force[i] += np.sum(fArr)
        return force


def calculateForceOnBoundary(blocks):
    force = np.array([0.0] * 3)
    for block in blocks:
        force += block['ForceCalculation'].calculateForceOnBoundary(block['pdfs'])
    return np.array(mpi.reduceReal(force, mpi.SUM))


def makeNacaAirfoilImage(domainSize, thicknessInPercent=30, angle=0):
    def nacaAirfoil(x, thicknessInPercent, chordLength):
        xOverC = x / chordLength
        y_t = 0
        coeffs = [0.2969, -0.1260, - 0.3516, 0.2843, -0.1015]
        for coeff, exponent in zip(coeffs, [0.5, 1, 2, 3, 4]):
            y_t += coeff * xOverC ** exponent
        y_t *= 5 * thicknessInPercent / 100 * chordLength
        return y_t

    domain = np.zeros(domainSize)
    it = np.nditer(domain, flags=['multi_index'], op_flags=['readwrite'])
    while not it.finished:
        x, y = it.multi_index
        y -= domain.shape[1] / 2
        if abs(y) < nacaAirfoil(x, thicknessInPercent, domain.shape[0]):
            it[0] = 1
        it.iternext()
    domain = np.rot90(domain, 1)
    domain = scipy.ndimage.interpolation.rotate(domain, angle=angle)

    domain[domain > 0.5] = 1
    domain[domain <= 0.5] = 0
    domain = domain.astype(np.int32)
    return domain


img = makeNacaAirfoilImage([300, 300], 30, angle=-30)

omega = 1.9
blocks = createUniformBlockGrid(cells=(500, 200, 1), periodic=(1, 0, 1))

collisionModel = lbm.collisionModels.SRT(omega)
forceModel = lbm.forceModels.SimpleConstant((1e-5, 0, 0))
latticeModel = lbm.makeLatticeModel("D2Q9", collisionModel, forceModel)
lbm.addPdfFieldToStorage(blocks, "pdfs", latticeModel, velocityAdaptor="vel", densityAdaptor="rho", initialDensity=1.0)
field.addFlagFieldToStorage(blocks, 'flags')
lbm.addBoundaryHandlingToStorage(blocks, 'boundary', 'pdfs', 'flags')

# setBoundaryFromArray( blocks, 'boundary', makeSlice[0.4:0.6, 0.4:0.55 ,0.5], img, { 1: 'NoSlip' } )
setBoundaryFromBlackAndWhiteImage(blocks, "boundary", makeSlice[0.25:0.75, 0.3:0.6, 0.5], imageFile, "NoSlip")
setBoundariesChannel(blocks, 'boundary')

blocks.addBlockData('ForceCalculation', ForceCalculationMasks.addToBlock)

sweep = lbm.makeCellwiseSweep(blocks, "pdfs", flagFieldID='flags', flagList=['fluid'])

scheme = createUniformBufferedScheme(blocks, 'D3Q19')
scheme.addDataToCommunicate(field.createPackInfo(blocks, 'pdfs'))


def timestep():
    scheme()
    for block in blocks:
        block['boundary']()
    for block in blocks:
        sweep.streamCollide(block)
    return calculateForceOnBoundary(blocks)


def run(timesteps):
    for t in range(timesteps):
        scheme()
        for block in blocks:
            block['boundary']()
        for block in blocks:
            sweep.streamCollide(block)


def makeAnimation(blocks, interval=30, frames=180):
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    plt.style.use('ggplot')
    NR_OF_TIMESTEPS_SHOWN = 600
    lifts = []

    fig = plt.gcf()
    f = field.gather(blocks, 'rho', makeSlice[:, :, 0.5])
    im = None

    ymax = [0.05]
    if f:
        npField = field.toArray(f).squeeze()
        npField = np.rot90(npField, 1)

        plt.subplot(2, 1, 1)
        plt.title("Lattice Density")
        im = plt.imshow(npField)
        plt.colorbar()

        plt.subplot(2, 1, 2)
        plt.title("Lift")
        plt.ylim(0, ymax[0])
        plt.xlim(0, NR_OF_TIMESTEPS_SHOWN)
        liftPlot, = plt.plot(lifts)

    def updatefig(*args):
        force = timestep()
        f = field.gather(blocks, 'rho', makeSlice[:, :, 0.5])
        if f:
            npField = field.toArray(f).squeeze()
            npField = np.rot90(npField, 1)
            im.set_array(npField)
            im.autoscale()
            if lifts and max(lifts) * 1.2 > ymax[0]:
                ymax[0] = max(lifts) * 1.2
                liftPlot.axes.set_ylim(0, ymax[0])

            lifts.append(force[1])
            nrOfSamples = len(lifts)
            xMin = max(0, nrOfSamples - NR_OF_TIMESTEPS_SHOWN)
            liftPlot.axes.set_xlim(xMin, xMin + NR_OF_TIMESTEPS_SHOWN)
            liftPlot.set_data(np.arange(nrOfSamples), lifts)
            return im, liftPlot

    return animation.FuncAnimation(fig, updatefig, interval=interval, frames=frames, blit=False, repeat=False)


showPlots = False

if showPlots:
    import waLBerla.plot as wplt

    setFieldUsingFlagMask(blocks, 'pdfs', np.NaN, 'flags', ['NoSlip'])
    run(1)
    setFieldUsingFlagMask(blocks, 'pdfs', np.NaN, 'flags', ['NoSlip'])

    ani = makeAnimation(blocks, frames=6000, )
    wplt.show()
else:
    run(10)
