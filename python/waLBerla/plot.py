"""Small wrappers aroung matplotlib to gather and plot parts of a waLBerla domain"""

import matplotlib.animation as animation
import numpy as np

try:
    from . import walberla_cpp
except ImportError:
    import walberla_cpp

from matplotlib.pyplot import imshow, gcf, figure, plot, quiver


def field_show(numpy_field, **kwargs):
    numpy_field = np.rot90(numpy_field, 3)
    imshow(numpy_field, origin='lower', **kwargs)


def scalar_field(blocks, name, slice_definition, f_coordinate=0, target_rank=0, **kwargs):
    """Plots a 2D slice through the global domain as an image

    :param blocks:              the blockforest
    :param name:                Name of the block data to be plotted. Has to be a scalar field
    :param slice_definition:    a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param f_coordinate:        value of the forth field coordinate (f)
    :param target_rank:         rank that gathers and plots the data
    :param kwargs:              further keyword arguments are passed to matplotlib.pyplot.imshow
    """
    f = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
    if f:
        numpy_field = np.asarray(f).squeeze()
        numpy_field = np.swapaxes(numpy_field, 0, 1)
        imshow(numpy_field, origin='lower', **kwargs)


def scalar_field_animation(blocks, name, slice_definition, run_function, plot_setup_function=lambda: None,
                           plot_update_function=lambda: None, f_coordinate=0, target_rank=0,
                           interval=30, frames=180, **kwargs):
    """Creates animation of 2D slices through the global domain

    :param blocks:               the blockforest
    :param name:                 Name of the block data to be plotted. Has to be a scalar field
    :param slice_definition:     a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param run_function:         function without arguments which is run between frames (should move simulation forward)
    :param plot_setup_function:  function without arguments that is called after the plot was initially created.
                                 Can be used to configure plot (set title etc.)
    :param plot_update_function: function without arguments that is called when figure is updated
    :param f_coordinate:         value of the forth field coordinate (f)
    :param target_rank:          rank that gathers and plots the data
    :param interval:             passed to matplotlib.animation.FuncAnimation: milliseconds between two frames
    :param frames:               passed to :class:`matplotlib.animation.FuncAnimation` number of frames

    for other params see :func:`scalarField`
    """
    fig = gcf()
    f = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
    im = None
    if f:
        numpy_field = np.asarray(f).squeeze()
        numpy_field = np.swapaxes(numpy_field, 0, 1)
        im = imshow(numpy_field, origin='lower', **kwargs)
        plot_setup_function()

    def updatefig(*args):
        run_function()
        gathered_field = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
        if gathered_field:
            n_field = np.swapaxes(np.asarray(gathered_field), 0, 1)
            if n_field.shape == 4:
                n_field = n_field[:, :, :, f_coordinate].squeeze()
            else:
                n_field = n_field.squeeze()
            im.set_array(n_field)
            plot_update_function()
            return im,

    return animation.FuncAnimation(fig, updatefig, interval=interval, frames=frames)


def vector_field(blocks, name, slice_definition, x_component=0, y_component=1, target_rank=0,
                 x_step=1, y_step=1, **kwargs):
    """Plots a vector field slice using matplotlib quiver

    :param blocks:           the blockforest
    :param name:             Name of the block data to be plotted. Has to be a scalar field
    :param slice_definition: a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param x_component:      which component of the vector field (0,1 or 2)
                             to take as the horizontal value for the quiver arrows
    :param y_component:      which component of the vector field (0,1 or 2)
                             to take as the vertical value for the quiver arrows
    :param x_step:           take only every xStep's cell/arrow in x direction
    :param y_step:           take only every yStep's cell/arrow in y direction
    :param target_rank:      rank that gathers and plots the data
    :param kwargs:           further keyword arguments are passed to matplotlib.pyplot.quiver
    """
    f = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
    if f:
        numpy_field = np.swapaxes(np.asarray(f), 0, 1)
        x_vel = numpy_field[::x_step, ::y_step, :, x_component].squeeze()
        y_vel = numpy_field[::x_step, ::y_step, :, y_component].squeeze()
        quiver(x_vel, y_vel, **kwargs)


def along_line(blocks, name, slice_definition, f_coordinate=0, target_rank=0, **kwargs):
    """Plot a field value along a one dimensional slice through the domain

    :param blocks:           the blockstorage
    :param name:             Name of the block data to be plotted. Has to be a scalar field
    :param slice_definition: a one dimensional slice through the domain. Can be created with :func:`waLBerla.makeSlice`
    :param f_coordinate:     value of the forth field coordinate (f)
    :param target_rank:      rank that gathers and plots the data
    :param kwargs:           further keyword arguments are passed to :func:`matplotlib.pyplot.plot`
    """
    f = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
    if f:
        npField = np.asarray(f)
        npField = npField[:, :, :, f_coordinate].squeeze()
        plot(npField, **kwargs)


def along_line_animation(blocks, name, slice_definition, run_function, plot_setup_function=lambda: None,
                         f_coordinate=0, target_rank=0, interval=30, frames=180, **kwargs):
    """Animated version of :func:`alongLine`

    :param blocks:              the blockforest
    :param name:                Name of the block data to be plotted. Has to be a scalar field
    :param slice_definition:    a two dimensional slice through the domain. Can be created with waLBerla.makeSlice
    :param run_function:        function without arguments which is run between frames (should move simulation forward)
    :param plot_setup_function: function without arguments that is called after the plot was initially created.
                                Can be used to configure plot (set title etc.)
    :param f_coordinate:        value of the forth field coordinate (f)
    :param target_rank:         rank that gathers and plots the data
    :param interval:            passed to matplotlib.animation.FuncAnimation: milliseconds between two frames
    :param frames:              passed to :class:`matplotlib.animation.FuncAnimation` number of frames

    """
    fig = figure()

    f = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
    line = None
    if f:
        numpy_field = np.asarray(f)
        numpy_field = numpy_field[:, :, :, f_coordinate].squeeze()
        line, = plot(numpy_field, **kwargs)
        plot_setup_function()

    def updatefig(*args):
        run_function()
        gathered_array = walberla_cpp.field.gather(blocks, name, slice_definition, targetRank=target_rank)
        if gathered_array:
            n_field = np.asarray(gathered_array)
            n_field = n_field[:, :, :, f_coordinate].squeeze()
            fig.gca().set_ylim((np.min(n_field), np.max(n_field)))
            line.set_ydata(n_field)
            return line,

    return animation.FuncAnimation(fig, updatefig, interval=interval, frames=frames, blit=False)
