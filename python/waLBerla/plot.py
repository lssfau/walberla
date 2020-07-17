"""Small wrappers aroung matplotlib to gather and plot parts of a waLBerla domain"""

import matplotlib.animation as animation
import numpy as np

try:
    from . import walberla_cpp
except ImportError:
    import walberla_cpp

from matplotlib.pyplot import imshow, gcf, figure, plot, quiver


def fieldShow(npField, **kwargs):
    npField = np.rot90(npField, 3)
    imshow(npField, origin='lower', **kwargs)


def scalarField(blocks, name, sliceDef, fCoord=0, targetRank=0, **kwargs):
    """Plots a 2D slice through the global domain as an image

    :param blocks:      the blockstorage
    :param name:        Name of the block data to be plotted. Has to be a scalar field
    :param sliceDef:    a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param fCoord:      value of the forth field coordinate (f)
    :param targetRank:  rank that gathers and plots the data
    :param kwargs:      further keyword arguments are passed to matplotlib.pyplot.imshow
    """
    f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
    if f:
        npField = np.asarray(f.buffer())[:, :, :, fCoord].squeeze()
        npField = np.swapaxes(npField, 0, 1)
        imshow(npField, origin='lower', **kwargs)


def scalarFieldAnimation(blocks, name, sliceDef, runFunction, plotSetupFunction=lambda: None,
                         plotUpdateFunction=lambda: None, fCoord=0, targetRank=0, interval=30, frames=180, **kwargs):
    """Creates animation of 2D slices through the global domain

    :param runFunction:        function without arguments which is run between frames (should move simulation forward)
    :param plotSetupFunction:  function without arguments that is called after the plot was initially created.
                               Can be used to configure plot (set title etc.)
    :param plotUpdateFunction: function without arguments that is called when figure is updated
    :param interval:           passed to matplotlib.animation.FuncAnimation: milliseconds between two frames
    :param frames:             passed to :class:`matplotlib.animation.FuncAnimation` number of frames

    for other params see :func:`scalarField`
    """
    fig = gcf()
    f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
    im = None
    if f:
        npField = np.asarray(f.buffer())[:, :, :, fCoord].squeeze()
        npField = np.swapaxes(npField, 0, 1)
        im = imshow(npField, origin='lower', **kwargs)
        plotSetupFunction()

    def updatefig(*args):
        runFunction()
        f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
        if f:
            npField = np.swapaxes(np.asarray(f.buffer()), 0, 1)
            npField = npField[:, :, :, fCoord].squeeze()
            im.set_array(npField)
            plotUpdateFunction()
            return im,

    return animation.FuncAnimation(fig, updatefig, interval=interval, frames=frames)


def vectorField(blocks, name, sliceDef, xComponent=0, yComponent=1, targetRank=0, xStep=1, yStep=1, **kwargs):
    """Plots a vector field slice using matplotlib quiver

    :param blocks:     the blockstorage
    :param name:       Name of the block data to be plotted. Has to be a scalar field
    :param sliceDef:   a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param xComponent: which component of the vector field (0,1 or 2)
                       to take as the horizontal value for the quiver arrows
    :param yComponent: which component of the vector field (0,1 or 2)
                       to take as the vertical value for the quiver arrows
    :param xStep:      take only every xStep's cell/arrow in x direction
    :param yStep:      take only every yStep's cell/arrow in y direction
    :param targetRank: rank that gathers and plots the data
    :param kwargs:     further keyword arguments are passed to matplotlib.pyplot.quiver
    """
    f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
    if f:
        npField = np.swapaxes(np.asarray(f.buffer()), 0, 1)
        xVel = npField[::xStep, ::yStep, :, xComponent].squeeze()
        yVel = npField[::xStep, ::yStep, :, yComponent].squeeze()
        quiver(xVel, yVel, **kwargs)


def alongLine(blocks, name, sliceDef, fCoord=0, targetRank=0, **kwargs):
    """Plot a field value along a one dimensional slice through the domain

    :param blocks:      the blockstorage
    :param name:        Name of the block data to be plotted. Has to be a scalar field
    :param sliceDef:    a one dimensional slice through the domain. Can be created with :func:`waLBerla.makeSlice`
    :param fCoord:      value of the forth field coordinate (f)
    :param targetRank:  rank that gathers and plots the data
    :param kwargs:      further keyword arguments are passed to :func:`matplotlib.pyplot.plot`
    """
    f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
    if f:
        npField = np.asarray(f.buffer())
        npField = npField[:, :, :, fCoord].squeeze()
        plot(npField, **kwargs)


def alongLineAnimation(blocks, name, sliceDef, runFunction, plotSetupFunction=lambda: None, fCoord=0, targetRank=0,
                       interval=30, frames=180, **kwargs):
    """Animated version of :func:`alongLine`

    For parameter documentation see :func:`scalarFieldAnimation` and :func:`alongLine`

    """
    fig = figure()

    f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
    line = None
    if f:
        npField = np.asarray(f.buffer())
        npField = npField[:, :, :, fCoord].squeeze()
        line, = plot(npField, **kwargs)
        plotSetupFunction()

    def updatefig(*args):
        runFunction()
        f = walberla_cpp.field.gather(blocks, name, sliceDef, targetRank=targetRank)
        if f:
            npField = np.asarray(f.buffer())
            npField = npField[:, :, :, fCoord].squeeze()
            fig.gca().set_ylim((np.min(npField), np.max(npField)))
            line.set_ydata(npField)
            return line,

    return animation.FuncAnimation(fig, updatefig, interval=interval, frames=frames, blit=False)
