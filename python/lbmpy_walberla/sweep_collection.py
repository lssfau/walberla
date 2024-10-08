from dataclasses import replace
from typing import Dict

from jinja2 import Environment, PackageLoader, StrictUndefined

import sympy as sp
import numpy as np

from pystencils import Target, create_kernel, Assignment
from pystencils.bit_masks import flag_cond
from pystencils.config import CreateKernelConfig
from pystencils.field import Field, fields
from pystencils.simp import add_subexpressions_for_field_reads
from pystencils.typing import BasicType, PointerType, FieldPointerSymbol, TypedSymbol, CastFunc

from lbmpy.advanced_streaming import is_inplace, get_accessor, Timestep
from lbmpy.creationfunctions import LbmCollisionRule, LBMConfig, LBMOptimisation
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter, macroscopic_values_getter
from lbmpy.updatekernels import create_lbm_kernel, create_stream_only_kernel

from pystencils_walberla.kernel_selection import KernelCallNode, KernelFamily
from pystencils_walberla.utility import config_from_context
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from lbmpy_walberla.utility import create_pdf_field, timestep_suffix

from .alternating_sweeps import EvenIntegerCondition
from .function_generator import kernel_family_function_generator


def generate_lbm_sweep_collection(ctx, class_name: str, collision_rule: LbmCollisionRule,
                                  lbm_config: LBMConfig, lbm_optimisation: LBMOptimisation,
                                  refinement_scaling=None, macroscopic_fields: Dict[str, Field] = None,
                                  target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None,
                                  max_threads=None, set_pre_collision_pdfs=True,
                                  **create_kernel_params):

    config = config_from_context(ctx, target=target, data_type=data_type,
                                 cpu_openmp=cpu_openmp, cpu_vectorize_info=cpu_vectorize_info, **create_kernel_params)

    streaming_pattern = lbm_config.streaming_pattern
    field_layout = lbm_optimisation.field_layout

    # usually a numpy layout is chosen by default i.e. xyzf - which is bad for waLBerla where at least the spatial
    # coordinates should be ordered in reverse direction i.e. zyx
    lb_method = collision_rule.method

    if field_layout == 'fzyx':
        config.cpu_vectorize_info['assume_inner_stride_one'] = True
    elif field_layout == 'zyxf':
        config.cpu_vectorize_info['assume_inner_stride_one'] = False

    src_field = lbm_optimisation.symbolic_field
    if not src_field:
        src_field = create_pdf_field(config=config, name="pdfs", stencil=lbm_config.stencil,
                                     field_layout=lbm_optimisation.field_layout)
    if is_inplace(streaming_pattern):
        dst_field = src_field
    else:
        dst_field = lbm_optimisation.symbolic_temporary_field
        if not dst_field:
            dst_field = create_pdf_field(config=config, name="pdfs_tmp", stencil=lbm_config.stencil,
                                         field_layout=lbm_optimisation.field_layout)

    config = replace(config, ghost_layers=0)
    function_generators = []

    def family(name):
        return lbm_kernel_family(class_name, name, collision_rule, streaming_pattern, src_field, dst_field, config)

    def generator(name, kernel_family):
        return kernel_family_function_generator(name, kernel_family, namespace='lbm', max_threads=max_threads)

    all_fields = collision_rule.bound_fields.union(collision_rule.free_fields)
    all_fields.update({src_field, dst_field})
    all_fields = list(sorted(all_fields, key=lambda e: str(e)))

    bw_stream_collide = block_wise_stream_collide(class_name, collision_rule, lbm_config, src_field, dst_field, config)
    bw_stream = block_wise_stream(class_name, collision_rule, lbm_config, src_field, dst_field, config)

    function_generators.append(generator('streamCollide', family("streamCollide")))
    function_generators.append(generator('collide', family("collide")))
    function_generators.append(generator('stream', family("stream")))
    function_generators.append(generator('streamOnlyNoAdvancement', family("streamOnlyNoAdvancement")))

    config_unoptimized = replace(config, cpu_vectorize_info=None, cpu_prepend_optimizations=[], cpu_blocking=None)

    setter_family = get_setter_family(class_name, lb_method, src_field, streaming_pattern, macroscopic_fields,
                                      config_unoptimized, set_pre_collision_pdfs)
    setter_generator = kernel_family_function_generator('initialise', setter_family,
                                                        namespace='lbm', max_threads=max_threads)
    function_generators.append(setter_generator)

    getter_family = get_getter_family(class_name, lb_method, src_field, streaming_pattern, macroscopic_fields,
                                      config_unoptimized)
    getter_generator = kernel_family_function_generator('calculateMacroscopicParameters', getter_family,
                                                        namespace='lbm', max_threads=max_threads)
    function_generators.append(getter_generator)

    contexts_function_generators = list()
    for fct in function_generators:
        contexts_function_generators.append(fct())

    namespaces = set([context['namespace'] for context in contexts_function_generators])
    assert len(namespaces) == 1, "All function_generators must output the same namespace!"
    namespace = namespaces.pop()

    headers = set()
    for context in contexts_function_generators:
        for header in context['interface_spec'].headers:
            headers.add(header)
        for header in context['kernel'].get_headers():
            headers.add(header)

    kernel_list = list()
    for context in contexts_function_generators:
        kernel_list.append(context['kernel'])

    kernels = list()
    for context in contexts_function_generators:
        kernels.append({
            'kernel': context['kernel'],
            'function_name': context['function_name'],
            'ghost_layers_to_include': 'ghost_layers',
            'field': context['field'],
            'max_threads': context['max_threads']
        })

    target = kernels[0]['kernel'].target

    jinja_context = {
        'block_stream_collide': bw_stream_collide,
        'block_stream': bw_stream,
        'all_fields': all_fields,
        'pdf_field': src_field,
        'kernel_list': kernel_list,
        'kernels': kernels,
        'namespace': namespace,
        'class_name': class_name,
        'headers': headers,
        'target': target.name.lower(),
        'is_gpu': target == Target.GPU,
        'parameter_scaling': refinement_scaling,
        'stencil_name': lbm_config.stencil.name,
        'D': lbm_config.stencil.D,
        'Q': lbm_config.stencil.Q,
        'inplace': is_inplace(lbm_config.streaming_pattern)
    }

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template("LBMSweepCollection.tmpl.h").render(**jinja_context)
    source = env.get_template("LBMSweepCollection.tmpl.cpp").render(**jinja_context)

    source_extension = "cu" if target == Target.GPU and ctx.cuda else "cpp"
    ctx.write_file(f"{class_name}.h", header)
    ctx.write_file(f"{class_name}.{source_extension}", source)


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

    default_dtype = config.data_type.default_factory()
    if kernel_name == "streamCollide":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_lbm_kernel(collision_rule, src_field, dst_field, field_accessor, data_type=default_dtype)
        advance_timestep = {"field_name": src_field.name, "function": "advanceTimestep"}
        temporary_fields = ['pdfs_tmp']
        field_swaps = [('pdfs', 'pdfs_tmp')]
    elif kernel_name == "collide":
        def lbm_kernel(field_accessor, lb_stencil):
            return create_lbm_kernel(collision_rule, src_field, dst_field, CollideOnlyInplaceAccessor(),
                                     data_type=default_dtype)
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


def get_setter_family(class_name, lb_method, pdfs, streaming_pattern, macroscopic_fields,
                      config: CreateKernelConfig, set_pre_collision_pdfs: bool):
    dim = lb_method.stencil.D
    density = macroscopic_fields.get('density', 1.0)
    velocity = macroscopic_fields.get('velocity', [0.0] * dim)

    default_dtype = config.data_type.default_factory()

    get_timestep = {"field_name": pdfs.name, "function": "getTimestepPlusOne"}
    temporary_fields = ()
    field_swaps = ()

    if is_inplace(streaming_pattern):
        nodes = list()
        for timestep in [Timestep.EVEN, Timestep.ODD]:
            timestep_suffix = str(timestep)
            setter = macroscopic_values_setter(lb_method,
                                               density=density, velocity=velocity, pdfs=pdfs,
                                               streaming_pattern=streaming_pattern, previous_timestep=timestep,
                                               set_pre_collision_pdfs=set_pre_collision_pdfs)

            if default_dtype != pdfs.dtype:
                setter = add_subexpressions_for_field_reads(setter, data_type=default_dtype)

            setter_ast = create_kernel(setter, config=config)
            setter_ast.function_name = 'kernel_initialise' + timestep_suffix
            nodes.append(KernelCallNode(setter_ast))
        tree = EvenIntegerCondition('timestep', nodes[0], nodes[1], parameter_dtype=np.uint8)
        family = KernelFamily(tree, class_name, field_timestep=get_timestep)
    else:
        timestep = Timestep.BOTH
        setter = macroscopic_values_setter(lb_method,
                                           density=density, velocity=velocity, pdfs=pdfs,
                                           streaming_pattern=streaming_pattern, previous_timestep=timestep,
                                           set_pre_collision_pdfs=set_pre_collision_pdfs)

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

    default_dtype = config.data_type.default_factory()

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

            if default_dtype != pdfs.dtype:
                getter = add_subexpressions_for_field_reads(getter, data_type=default_dtype)

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


def block_wise_stream_collide(class_name, collision_rule, lbm_config, src_field, dst_field, config):

    if not is_inplace(lbm_config.streaming_pattern):
        return None
    else:
        ast_even, all_fields = create_block_wise_ast(collision_rule, src_field, dst_field,
                                                     lbm_config, Timestep.EVEN, config, False)
        even_call = KernelCallNode(ast_even)
        ast_odd, _ = create_block_wise_ast(collision_rule, src_field, dst_field,
                                           lbm_config, Timestep.ODD, config, False)
        odd_call = KernelCallNode(ast_odd)
        tree = EvenIntegerCondition('timestep', even_call, odd_call, parameter_dtype=np.uint8)

    family = KernelFamily(tree, class_name)

    indexed_to_field_name = dict()
    for field in all_fields:
        indexed_to_field_name[field.name] = f"_data_{field.name}_dp"

    context = {
        'kernel': family,
        'all_fields': all_fields,
        'namespace': 'lbm',
        'function_name': 'blockStreamCollide',
        'indexed_to_field_name': indexed_to_field_name,
        'max_threads': None
    }

    return context


def block_wise_stream(class_name, collision_rule, lbm_config, src_field, dst_field, config):

    if not is_inplace(lbm_config.streaming_pattern):
        return None
    else:
        ast_even, all_fields = create_block_wise_ast(collision_rule, src_field, dst_field,
                                                     lbm_config, Timestep.EVEN, config, True)
        even_call = KernelCallNode(ast_even)
        ast_odd, _ = create_block_wise_ast(collision_rule, src_field, dst_field,
                                           lbm_config, Timestep.ODD, config, True)
        odd_call = KernelCallNode(ast_odd)
        tree = EvenIntegerCondition('timestep', even_call, odd_call, parameter_dtype=np.uint8)

    family = KernelFamily(tree, class_name)

    indexed_to_field_name = dict()
    for field in all_fields:
        indexed_to_field_name[field.name] = f"_data_{field.name}_dp"

    context = {
        'kernel': family,
        'all_fields': all_fields,
        'namespace': 'lbm',
        'function_name': 'blockStream',
        'indexed_to_field_name': indexed_to_field_name,
        'max_threads': None
    }

    return context


def create_block_wise_ast(collision_rule, src_field, dst_field, lbm_config, timestep, config, stream_only):
    stencil = lbm_config.stencil
    streaming_pattern = lbm_config.streaming_pattern
    default_dtype = config.data_type.default_factory()
    config = replace(config, gpu_indexing_params={})

    accessor = get_accessor(streaming_pattern, timestep)

    if stream_only:
        update_rule = create_stream_only_kernel(stencil, src_field, dst_field, accessor)
    else:
        update_rule = create_lbm_kernel(collision_rule, src_field, dst_field, accessor, data_type=default_dtype)

    bound_fields = update_rule.bound_fields
    free_fields = update_rule.free_fields

    all_fields = list(bound_fields.union(free_fields))
    all_fields.sort(key=lambda field: field.name)

    index = TypedSymbol("index", dtype=BasicType(np.int64))
    index_shape = TypedSymbol("_size_0", dtype=BasicType(np.int64))

    ass = list()
    for field in all_fields:
        const = True if field in free_fields else False
        ptr_type = PointerType(field.dtype, const=const, restrict=True, double_pointer=True)
        ptr = FieldPointerSymbol(field.name, field.dtype, const=const)
        f = sp.IndexedBase(TypedSymbol(f"_data_{field.name}_dp", dtype=ptr_type), shape=index_shape)
        ass.append(Assignment(ptr, f[index]))

    update_rule = ass + update_rule.all_assignments

    ast = create_kernel(update_rule, config=config)
    base_name = "kernel_BlockStream" if stream_only else "kernel_BlockStreamCollide"
    ast.function_name = base_name + timestep_suffix(timestep)
    ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']
    return ast, all_fields
