import numpy as np
from pystencils import CreateKernelConfig, fields

from lbmpy.advanced_streaming import Timestep
from lbmpy.stencils import LBStencil


def timestep_suffix(timestep: Timestep):
    """ get the suffix as string for a timestep

    :param timestep: instance of class lbmpy.advanced_streaming.Timestep
    :return: either "even", "odd" or an empty string
    """
    return ("_" + str(timestep)) if timestep != Timestep.BOTH else ''


def create_pdf_field(config: CreateKernelConfig, name: str, stencil: LBStencil, field_layout: str = 'fzyx'):
    default_dtype = config.data_type.default_factory()
    data_type = default_dtype.numpy_dtype
    return fields(f'{name}({stencil.Q}) :{data_type}[{stencil.D}D]', layout=field_layout)

