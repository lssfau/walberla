from dataclasses import replace
from typing import Dict

import sympy as sp
import numpy as np

from pystencils import Target, create_kernel
from pystencils.config import CreateKernelConfig
from pystencils.field import Field

from lbmpy.advanced_streaming import is_inplace, get_accessor, Timestep
from lbmpy.creationfunctions import LbmCollisionRule
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter, macroscopic_values_getter
from lbmpy.updatekernels import create_lbm_kernel, create_stream_only_kernel

from pystencils_walberla.kernel_selection import KernelCallNode, KernelFamily
from pystencils_walberla.utility import config_from_context
from pystencils_walberla import generate_sweep_collection

from .alternating_sweeps import EvenIntegerCondition
from .function_generator import kernel_family_function_generator


def generate_lbm_sweep_collection(ctx, class_name: str, collision_rule: LbmCollisionRule,
                                  streaming_pattern='pull',
                                  field_layout='fzyx', refinement_scaling=None,
                                  macroscopic_fields: Dict[str, Field] = None,
                                  target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None,
                                  max_threads=None,
                                  **create_kernel_params):
    config = config_from_context(ctx, target=target, data_type=data_type,
                                 cpu_openmp=cpu_openmp, cpu_vectorize_info=cpu_vectorize_info, **create_kernel_params)

    # usually a numpy layout is chosen by default i.e. xyzf - which is bad for waLBerla where at least the spatial
    # coordinates should be ordered in reverse direction i.e. zyx
    lb_method = collision_rule.method

    q = lb_method.stencil.Q
    dim = lb_method.stencil.D

    if field_layout == 'fzyx':
        config.cpu_vectorize_info['assume_inner_stride_one'] = True
    elif field_layout == 'zyxf':
        config.cpu_vectorize_info['assume_inner_stride_one'] = False

    src_field = Field.create_generic('pdfs', dim, config.data_type['pdfs'].numpy_dtype,
                                     index_dimensions=1, layout=field_layout, index_shape=(q,))
    if is_inplace(streaming_pattern):
        dst_field = src_field
    else:
        dst_field = Field.create_generic('pdfs_tmp', dim, config.data_type['pdfs_tmp'].numpy_dtype,
                                         index_dimensions=1, layout=field_layout,
                                         index_shape=(q,))

    config = replace(config, ghost_layers=0)

    function_generators = []

    def family(name):
        return lbm_kernel_family(class_name, name, collision_rule, streaming_pattern, src_field, dst_field, config)

    def generator(name, kernel_family):
        return kernel_family_function_generator(name, kernel_family, namespace='lbm', max_threads=max_threads)

    function_generators.append(generator('streamCollide', family("streamCollide")))
    function_generators.append(generator('collide', family("collide")))
    function_generators.append(generator('stream', family("stream")))
    function_generators.append(generator('streamOnlyNoAdvancement', family("streamOnlyNoAdvancement")))

    setter_family = get_setter_family(class_name, lb_method, src_field, streaming_pattern, macroscopic_fields, config)
    setter_generator = kernel_family_function_generator('initialise', setter_family,
                                                        namespace='lbm', max_threads=max_threads)
    function_generators.append(setter_generator)

    getter_family = get_getter_family(class_name, lb_method, src_field, streaming_pattern, macroscopic_fields, config)
    getter_generator = kernel_family_function_generator('calculateMacroscopicParameters', getter_family,
                                                        namespace='lbm', max_threads=max_threads)
    function_generators.append(getter_generator)

    generate_sweep_collection(ctx, class_name, function_generators, refinement_scaling)


class RefinementScaling:
    def __init__(self):
        self.scaling_info = []

    def add_standard_relaxation_rate_scaling(self, viscosity_relaxation_rate):
        self.add_scaling(viscosity_relaxation_rate)

    def add_scaling(self, parameter):
        if isinstance(parameter, sp.Symbol):
            self.scaling_info.append(parameter.name)
        else:
            raise ValueError("Only pure symbols allowed")


def lbm_kernel_family(class_name, kernel_name,
                      collision_rule, streaming_pattern, src_field, dst_field, config: CreateKernelConfig):

    if kernel_name == "streamCollide":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_lbm_kernel(collision_rule, src_field, dst_field, field_accessor)
        advance_timestep = {"field_name": src_field.name, "function": "advanceTimestep"}
        temporary_fields = ['pdfs_tmp']
        field_swaps = [('pdfs', 'pdfs_tmp')]
    elif kernel_name == "collide":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_lbm_kernel(collision_rule, src_field, dst_field, CollideOnlyInplaceAccessor())
        advance_timestep = {"field_name": src_field.name, "function": "advanceTimestep"}
        temporary_fields = ()
        field_swaps = ()
    elif kernel_name == "stream":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_stream_only_kernel(lb_stencil, src_field, dst_field, field_accessor)
        advance_timestep = {"field_name": src_field.name, "function": "advanceTimestep"}
        temporary_fields = ['pdfs_tmp']
        field_swaps = [('pdfs', 'pdfs_tmp')]
    elif kernel_name == "streamOnlyNoAdvancement":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_stream_only_kernel(lb_stencil, src_field, dst_field, field_accessor)
        advance_timestep = {"field_name": src_field.name, "function": "getTimestepPlusOne"}
        temporary_fields = ['pdfs_tmp']
        field_swaps = ()
    else:
        raise ValueError(f"kernel name: {kernel_name} is not valid")

    lb_method = collision_rule.method
    stencil = lb_method.stencil

    if is_inplace(streaming_pattern):
        nodes = list()
        for timestep in [Timestep.EVEN, Timestep.ODD]:
            accessor = get_accessor(streaming_pattern, timestep)
            timestep_suffix = str(timestep)

            update_rule = lbm_kernel(accessor, stencil)
            ast = create_kernel(update_rule, config=config)
            ast.function_name = 'kernel_' + kernel_name + timestep_suffix
            ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']
            nodes.append(KernelCallNode(ast))

        tree = EvenIntegerCondition('timestep', nodes[0], nodes[1], parameter_dtype=np.uint8)
        family = KernelFamily(tree, class_name, field_timestep=advance_timestep)
    else:
        timestep = Timestep.BOTH
        accessor = get_accessor(streaming_pattern, timestep)

        update_rule = lbm_kernel(accessor, stencil)
        ast = create_kernel(update_rule, config=config)
        ast.function_name = 'kernel_' + kernel_name
        ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']
        node = KernelCallNode(ast)
        family = KernelFamily(node, class_name, temporary_fields=temporary_fields, field_swaps=field_swaps)

    return family


def get_setter_family(class_name, lb_method, pdfs, streaming_pattern, macroscopic_fields, config: CreateKernelConfig):
    dim = lb_method.stencil.D
    density = macroscopic_fields.get('density', 1.0)
    velocity = macroscopic_fields.get('velocity', [0.0] * dim)

    get_timestep = {"field_name": pdfs.name, "function": "getTimestep"}
    temporary_fields = ()
    field_swaps = ()

    if is_inplace(streaming_pattern):
        nodes = list()
        for timestep in [Timestep.EVEN, Timestep.ODD]:
            timestep_suffix = str(timestep)
            setter = macroscopic_values_setter(lb_method,
                                               density=density, velocity=velocity, pdfs=pdfs,
                                               streaming_pattern=streaming_pattern, previous_timestep=timestep)

            setter_ast = create_kernel(setter, config=config)
            setter_ast.function_name = 'kernel_initialise' + timestep_suffix
            nodes.append(KernelCallNode(setter_ast))
        tree = EvenIntegerCondition('timestep', nodes[0], nodes[1], parameter_dtype=np.uint8)
        family = KernelFamily(tree, class_name, field_timestep=get_timestep)
    else:
        timestep = Timestep.BOTH
        setter = macroscopic_values_setter(lb_method,
                                           density=density, velocity=velocity, pdfs=pdfs,
                                           streaming_pattern=streaming_pattern, previous_timestep=timestep)

        setter_ast = create_kernel(setter, config=config)
        setter_ast.function_name = 'kernel_initialise'
        node = KernelCallNode(setter_ast)
        family = KernelFamily(node, class_name, temporary_fields=temporary_fields, field_swaps=field_swaps)

    return family


def get_getter_family(class_name, lb_method, pdfs, streaming_pattern, macroscopic_fields, config: CreateKernelConfig):
    density = macroscopic_fields.get('density', None)
    velocity = macroscopic_fields.get('velocity', None)

    if density is None and velocity is None:
        return None

    get_timestep = {"field_name": pdfs.name, "function": "getTimestep"}
    temporary_fields = ()
    field_swaps = ()

    if is_inplace(streaming_pattern):
        nodes = list()
        for timestep in [Timestep.EVEN, Timestep.ODD]:
            timestep_suffix = str(timestep)
            getter = macroscopic_values_getter(lb_method,
                                               density=density, velocity=velocity, pdfs=pdfs,
                                               streaming_pattern=streaming_pattern, previous_timestep=timestep)

            getter_ast = create_kernel(getter, config=config)
            getter_ast.function_name = 'kernel_getter' + timestep_suffix
            nodes.append(KernelCallNode(getter_ast))
        tree = EvenIntegerCondition('timestep', nodes[0], nodes[1], parameter_dtype=np.uint8)
        family = KernelFamily(tree, class_name, field_timestep=get_timestep)
    else:
        timestep = Timestep.BOTH
        getter = macroscopic_values_getter(lb_method,
                                           density=density, velocity=velocity, pdfs=pdfs,
                                           streaming_pattern=streaming_pattern, previous_timestep=timestep)

        getter_ast = create_kernel(getter, config=config)
        getter_ast.function_name = 'kernel_getter'
        node = KernelCallNode(getter_ast)
        family = KernelFamily(node, class_name, temporary_fields=temporary_fields, field_swaps=field_swaps)

    return family
