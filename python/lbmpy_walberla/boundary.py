import pystencils_walberla.boundary
from lbmpy.boundaries.boundaryhandling import create_lattice_boltzmann_boundary_kernel
from lbmpy.advanced_streaming import Timestep, is_inplace

from pystencils_walberla.kernel_selection import KernelCallNode
from lbmpy_walberla.alternating_sweeps import EvenIntegerCondition, OddIntegerCondition, TimestepTrackerMapping
from lbmpy_walberla.additional_data_handler import default_additional_data_handler

from pystencils.data_types import TypedSymbol

import numpy as np


def generate_boundary(generation_context,
                      class_name,
                      boundary_object,
                      lb_method,
                      field_name='pdfs',
                      streaming_pattern='pull',
                      prev_timestep=Timestep.BOTH,
                      additional_data_handler=None,
                      namespace='lbm',
                      **create_kernel_params):
    if boundary_object.additional_data and additional_data_handler is None:
        target = create_kernel_params.get('target', 'cpu')
        additional_data_handler = default_additional_data_handler(boundary_object, lb_method, field_name, target=target)

    def boundary_creation_function(field, index_field, stencil, boundary_functor, target='cpu', **kwargs):
        return create_lattice_boltzmann_boundary_kernel(field, index_field, lb_method, boundary_functor,
                                                        streaming_pattern=streaming_pattern,
                                                        prev_timestep=prev_timestep,
                                                        target=target,
                                                        **kwargs)

    pystencils_walberla.boundary.generate_boundary(generation_context,
                                                   class_name,
                                                   boundary_object,
                                                   field_name=field_name,
                                                   neighbor_stencil=lb_method.stencil,
                                                   index_shape=[len(lb_method.stencil)],
                                                   kernel_creation_function=boundary_creation_function,
                                                   namespace=namespace,
                                                   additional_data_handler=additional_data_handler,
                                                   **create_kernel_params)


def generate_alternating_lbm_boundary(generation_context,
                                      class_name,
                                      boundary_object,
                                      lb_method,
                                      field_name='pdfs',
                                      streaming_pattern='pull',
                                      after_collision=True,
                                      additional_data_handler=None,
                                      namespace='lbm',
                                      **create_kernel_params):
    if boundary_object.additional_data and additional_data_handler is None:
        target = create_kernel_params.get('target', 'cpu')
        additional_data_handler = default_additional_data_handler(boundary_object, lb_method, field_name, target=target)

    timestep_param_name = 'timestep'
    timestep_param_dtype = np.uint8
    timestep_param = TypedSymbol(timestep_param_name, timestep_param_dtype)

    def boundary_creation_function(field, index_field, stencil, boundary_functor, target='cpu', **kwargs):
        pargs = (field, index_field, lb_method, boundary_functor)
        kwargs = {'target': target, **kwargs}
        ast_even = create_lattice_boltzmann_boundary_kernel(*pargs,
                                                            streaming_pattern=streaming_pattern,
                                                            prev_timestep=Timestep.EVEN,
                                                            **kwargs)
        ast_even.function_name = 'even'
        kernel_even = KernelCallNode(ast_even)

        if is_inplace(streaming_pattern):
            ast_odd = create_lattice_boltzmann_boundary_kernel(*pargs,
                                                               streaming_pattern=streaming_pattern,
                                                               prev_timestep=Timestep.ODD,
                                                               **kwargs)
            ast_odd.function_name = 'odd'
            kernel_odd = KernelCallNode(ast_odd)
        else:
            kernel_odd = kernel_even

        if after_collision:
            return EvenIntegerCondition(timestep_param_name, kernel_even, kernel_odd, timestep_param_dtype)
        else:
            return OddIntegerCondition(timestep_param_name, kernel_even, kernel_odd, timestep_param_dtype)

    interface_mappings = [TimestepTrackerMapping(timestep_param)]

    pystencils_walberla.boundary.generate_boundary(generation_context,
                                                   class_name,
                                                   boundary_object,
                                                   field_name=field_name,
                                                   neighbor_stencil=lb_method.stencil,
                                                   index_shape=[len(lb_method.stencil)],
                                                   kernel_creation_function=boundary_creation_function,
                                                   namespace=namespace,
                                                   additional_data_handler=additional_data_handler,
                                                   interface_mappings=interface_mappings,
                                                   **create_kernel_params)
