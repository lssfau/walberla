import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
from tempfile import NamedTemporaryFile
import base64

from IPython import get_ipython

ipython = get_ipython()

ipython.magic("matplotlib inline")  # Show plots as images embedded in iPython notebook


def setMplFigureSize():
    matplotlib.rcParams['figure.figsize'] = (15.0, 12.0)


VIDEO_TAG = """<video controls width="80%">
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""


def __anim_to_html(anim, fps):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=fps, extra_args=['-vcodec', 'libx264', '-pix_fmt',
                                                   'yuv420p', '-profile:v', 'baseline', '-level', '3.0'])
            video = open(f.name, "rb").read()
        anim._encoded_video = base64.b64encode(video).decode('ascii')

    return VIDEO_TAG.format(anim._encoded_video)


def makeImshowAnimation(grid, gridUpdateFunction, frames=90, **kwargs):
    from functools import partial
    fig = plt.figure()
    im = plt.imshow(grid, interpolation='none')

    def updatefig(*args, **kwargs):
        image = kwargs['image']
        image = gridUpdateFunction(image)
        im.set_array(image)
        return im,

    return animation.FuncAnimation(fig, partial(updatefig, image=grid), frames=frames)


# -------   Version 1: Embed the animation as HTML5 video --------- ----------------------------------

def displayAsHtmlVideo(anim, fps=30, show=True, **kwargs):
    try:
        plt.close(anim._fig)
        res = __anim_to_html(anim, fps)
        if show:
            return HTML(res)
        else:
            return HTML("")
    except KeyboardInterrupt:
        pass


# -------   Version 2: Animation is shown in extra matplotlib window ----------------------------------


def displayInExtraWindow(animation, *args, **kwargs):
    fig = plt.gcf()
    try:
        fig.canvas.manager.window.raise_()
    except Exception:
        pass
    plt.show()


# -------   Version 3: Animation is shown in images that are updated directly in website --------------

def displayAsHtmlImage(animation, show=True, iterations=10000, *args, **kwargs):
    from IPython import display

    try:
        if show:
            fig = plt.gcf()
        if show:
            animation._init_draw()
        for i in range(iterations):
            if show:
                display.display(fig)
            animation._step()
            if show:
                display.clear_output(wait=True)
    except KeyboardInterrupt:
        display.clear_output(wait=False)


# Dispatcher

animation_display_mode = 'imageupdate'
display_animation_func = None


def disp(*args, **kwargs):
    if not display_animation_func:
        raise ("Call set_display_mode first")
    return display_animation_func(*args, **kwargs)


def set_display_mode(mode):
    from IPython import get_ipython
    ipython = get_ipython()
    global animation_display_mode
    global display_animation_func
    animation_display_mode = mode
    if animation_display_mode == 'video':
        ipython.magic("matplotlib inline")
        display_animation_func = displayAsHtmlVideo
    elif animation_display_mode == 'window':
        ipython.magic("matplotlib qt")
        display_animation_func = displayInExtraWindow
    elif animation_display_mode == 'imageupdate':
        ipython.magic("matplotlib inline")
        display_animation_func = displayAsHtmlImage
    else:
        raise Exception("Unknown mode. Available modes 'imageupdate', 'video' and 'window' ")


set_display_mode('imageupdate')
setMplFigureSize()
