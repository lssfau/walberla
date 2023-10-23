from jinja2 import Environment, PackageLoader, StrictUndefined

import pystencils_walberla.boundary
from lbmpy.boundaries.boundaryconditions import LbBoundary
from lbmpy.boundaries.boundaryhandling import create_lattice_boltzmann_boundary_kernel
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from lbmpy.advanced_streaming import Timestep, is_inplace

from pystencils_walberla.kernel_selection import KernelCallNode
from lbmpy_walberla.alternating_sweeps import EvenIntegerCondition, OddIntegerCondition, TimestepTrackerMapping
from lbmpy_walberla.additional_data_handler import default_additional_data_handler

from pystencils import Target

import numpy as np


def lbm_boundary_generator(class_name: str, flag_uid: str, boundary_object: LbBoundary, additional_data_handler=None,
                           field_data_type=None):
    def generation_function(ctx, lb_method, field_name='pdfs', spatial_shape=None,
                            streaming_pattern='pull', after_collision=True,
                            namespace='lbm',
                            **create_kernel_params):
        context = __generate_alternating_lbm_boundary(generation_context=ctx,
                                                      class_name=class_name,
                                                      boundary_object=boundary_object,
                                                      lb_method=lb_method,
                                                      field_name=field_name,
                                                      spatial_shape=spatial_shape,
                                                      field_data_type=field_data_type,
                                                      streaming_pattern=streaming_pattern,
                                                      after_collision=after_collision,
                                                      additional_data_handler=additional_data_handler,
                                                      namespace=namespace,
                                                      **create_kernel_params)

        return context

    return {'flag_id': flag_uid, 'generator': generation_function}


def generate_boundary_collection(generation_context,
                                 class_name,
                                 boundary_generators,
                                 lb_method,
                                 field_name='pdfs',
                                 spatial_shape=None,
                                 streaming_pattern='pull',
                                 prev_timestep=Timestep.BOTH,
                                 namespace='lbm',
                                 **create_kernel_params):

    kernel_list = []
    includes = []
    boundary_classes = []
    additional_data_handlers = []
    flag_uids = []
    object_names = []
    targets = []

    for boundary_generator in boundary_generators:
        boundary_functor = boundary_generator['generator']
        context = boundary_functor(generation_context, lb_method, field_name, spatial_shape,
                                   streaming_pattern, prev_timestep, namespace, **create_kernel_params)

        kernel_list.append(context['kernel'])
        includes.append(f"\"{context['class_name']}.h\"")
        boundary_classes.append(f"{context['namespace']}::{context['class_name']}")
        additional_data_handlers.append(context['additional_data_handler'])
        flag_uids.append(boundary_generator['flag_id'])
        object_names.append(f"{context['class_name']}Object")
        targets.append(f"{context['target']}")

    additional_constructor_arguments = [a.constructor_arguments[2:] for a in additional_data_handlers]

    assert len(set(targets)) == 1
    target = targets[0]

    jinja_context = {
        'kernel_list': kernel_list,
        'class_name': class_name,
        'target': target,
        'namespace': namespace,
        'includes': includes,
        'boundary_classes': boundary_classes,
        'additional_data_handlers': additional_data_handlers,
        'additional_constructor_arguments': additional_constructor_arguments,
        'flag_uids': flag_uids,
        'object_names': object_names
    }

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    env.globals.update(zip=zip)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template("BoundaryCollection.tmpl.h").render(**jinja_context)

    generation_context.write_file(f"{class_name}.h", header)


# Internal
def __generate_alternating_lbm_boundary(generation_context,
                                        class_name,
                                        boundary_object,
                                        lb_method,
                                        field_name='pdfs',
                                        spatial_shape=None,
                                        field_data_type=None,
                                        streaming_pattern='pull',
                                        after_collision=True,
                                        additional_data_handler=None,
                                        namespace='lbm',
                                        **create_kernel_params):
    if boundary_object.additional_data and additional_data_handler is None:
        target = create_kernel_params.get('target', Target.CPU)
        additional_data_handler = default_additional_data_handler(boundary_object, lb_method, field_name, target=target)

    timestep_param_name = 'timestep'
    timestep_param_dtype = np.uint8

    def boundary_creation_function(field, index_field, stencil, boundary_functor, target=Target.CPU, **kwargs):
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

    timestep_advancement = {"field_name": field_name, "function": "getTimestep"}
    context = pystencils_walberla.boundary.generate_boundary(generation_context,
                                                             class_name=class_name,
                                                             boundary_object=boundary_object,
                                                             field_name=field_name,
                                                             neighbor_stencil=lb_method.stencil,
                                                             index_shape=[lb_method.stencil.Q],
                                                             spatial_shape=spatial_shape,
                                                             field_data_type=field_data_type,
                                                             kernel_creation_function=boundary_creation_function,
                                                             namespace=namespace,
                                                             additional_data_handler=additional_data_handler,
                                                             field_timestep=timestep_advancement,
                                                             **create_kernel_params)
    return context
