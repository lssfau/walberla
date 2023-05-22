from typing import Callable, Sequence

from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils import Target, Assignment
from pystencils import Field, create_kernel, create_staggered_kernel
from pystencils.astnodes import KernelFunction

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.kernel_selection import KernelCallNode, KernelFamily, HighLevelInterfaceSpec
from pystencils_walberla.utility import config_from_context


def generate_sweep(ctx: CodeGenerationContext, class_name: str, assignments: Sequence[Assignment],
                   namespace: str = 'pystencils', field_swaps=(), staggered=False, varying_parameters=(),
                   inner_outer_split=False, ghost_layers_to_include=0,
                   target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None, max_threads=None,
                   **create_kernel_params):
    """Generates a waLBerla sweep from a pystencils representation.

    The constructor of the C++ sweep class expects all kernel parameters (fields and parameters) in alphabetical order.
    Fields have to passed using BlockDataID's pointing to walberla fields

    Args:
        ctx: build system context filled with information from waLBerla's CMake. The context for example
                            defines where to write generated files, if OpenMP is available or which SIMD instruction
                            set should be used. See waLBerla examples on how to get a context.
        class_name: name of the generated sweep class
        assignments: list of assignments defining the stencil update rule or a :class:`KernelFunction`
        namespace: the generated class is accessible as walberla::<namespace>::<class_name>
        field_swaps: sequence of field pairs (field, temporary_field). The generated sweep only gets the first field
                     as argument, creating a temporary field internally which is swapped with the first field after
                     each iteration.
        staggered: set to True to create staggered kernels with `pystencils.create_staggered_kernel`
        varying_parameters: Depending on the configuration, the generated kernels may receive different arguments for
                            different setups. To not have to adapt the C++ application when then parameter change,
                            the varying_parameters sequence can contain parameter names, which are always expected by
                            the C++ class constructor even if the kernel does not need them.
        inner_outer_split: if True generate a sweep that supports separate iteration over inner and outer regions
                           to allow for communication hiding.
        ghost_layers_to_include: determines how many ghost layers should be included for the Sweep.
                                 This is relevant if a setter kernel should also set correct values to the ghost layers.
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
        cpu_vectorize_info: dictionary containing necessary information for the usage of a SIMD instruction set.
        max_threads: only relevant for GPU kernels. Will be argument of `__launch_bounds__`
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    if staggered:
        assert 'omp_single_loop' not in create_kernel_params
        create_kernel_params['omp_single_loop'] = False
    config = config_from_context(ctx, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 cpu_vectorize_info=cpu_vectorize_info, **create_kernel_params)

    if isinstance(assignments, KernelFunction):
        ast = assignments
        target = ast.target
    elif not staggered:
        ast = create_kernel(assignments, config=config)
    else:
        # This should not be necessary but create_staggered_kernel does not take a config at the moment ...
        ast = create_staggered_kernel(assignments, **config.__dict__)

    ast.function_name = class_name.lower()

    selection_tree = KernelCallNode(ast)
    generate_selective_sweep(ctx, class_name, selection_tree, target=target, namespace=namespace,
                             field_swaps=field_swaps, varying_parameters=varying_parameters,
                             inner_outer_split=inner_outer_split, ghost_layers_to_include=ghost_layers_to_include,
                             cpu_vectorize_info=config.cpu_vectorize_info,
                             cpu_openmp=config.cpu_openmp, max_threads=max_threads)


def generate_selective_sweep(ctx, class_name, selection_tree, interface_mappings=(), target=None,
                             namespace='pystencils', field_swaps=(), varying_parameters=(),
                             inner_outer_split=False, ghost_layers_to_include=0,
                             cpu_vectorize_info=None, cpu_openmp=False, max_threads=None):
    """Generates a selective sweep from a kernel selection tree. A kernel selection tree consolidates multiple
    pystencils ASTs in a tree-like structure. See also module `pystencils_walberla.kernel_selection`.

    Args:
        ctx: see documentation of `generate_sweep`
        class_name: name of the generated sweep class
        selection_tree: Instance of `AbstractKernelSelectionNode`, root of the selection tree
        interface_mappings: sequence of `AbstractInterfaceArgumentMapping` instances for selection arguments of
                            the selection tree
        target: `None`, `Target.CPU` or `Target.GPU`; inferred from kernels if `None` is given.
        namespace: see documentation of `generate_sweep`
        field_swaps: see documentation of `generate_sweep`
        varying_parameters: see documentation of `generate_sweep`
        inner_outer_split: see documentation of `generate_sweep`
        ghost_layers_to_include: see documentation of `generate_sweep`
        cpu_vectorize_info: Dictionary containing information about CPU vectorization applied to the kernels
        cpu_openmp: Whether or not CPU kernels use OpenMP parallelization
        max_threads: only relevant for GPU kernels. Will be argument of `__launch_bounds__`
    """
    def to_name(f):
        return f.name if isinstance(f, Field) else f

    field_swaps = tuple((to_name(e[0]), to_name(e[1])) for e in field_swaps)
    temporary_fields = tuple(e[1] for e in field_swaps)

    kernel_family = KernelFamily(selection_tree, class_name,
                                 temporary_fields, field_swaps, varying_parameters)

    if target is None:
        target = kernel_family.get_ast_attr('target')
    elif target != kernel_family.get_ast_attr('target'):
        raise ValueError('Mismatch between target parameter and AST targets.')

    if not ctx.gpu and target == Target.GPU:
        return

    representative_field = {p.field_name for p in kernel_family.parameters if p.is_field_parameter}
    representative_field = sorted(representative_field)[0]

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    interface_spec = HighLevelInterfaceSpec(kernel_family.kernel_selection_parameters, interface_mappings)

    jinja_context = {
        'kernel': kernel_family,
        'namespace': namespace,
        'class_name': class_name,
        'target': target.name.lower(),
        'field': representative_field,
        'ghost_layers_to_include': ghost_layers_to_include,
        'inner_outer_split': inner_outer_split,
        'interface_spec': interface_spec,
        'generate_functor': True,
        'cpu_vectorize_info': cpu_vectorize_info,
        'cpu_openmp': cpu_openmp,
        'max_threads': max_threads
    }
    header = env.get_template("Sweep.tmpl.h").render(**jinja_context)
    source = env.get_template("Sweep.tmpl.cpp").render(**jinja_context)

    source_extension = "cpp" if target == Target.CPU else "cu"
    ctx.write_file(f"{class_name}.h", header)
    ctx.write_file(f"{class_name}.{source_extension}", source)


def generate_sweep_collection(ctx, class_name: str, function_generators: Sequence[Callable], parameter_scaling=None):
    """Generates a sweep collection
    """

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
        'kernel_list': kernel_list,
        'kernels': kernels,
        'namespace': namespace,
        'class_name': class_name,
        'headers': headers,
        'target': target.name.lower(),
        'parameter_scaling': parameter_scaling,
    }

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template("SweepCollection.tmpl.h").render(**jinja_context)
    source = env.get_template("SweepCollection.tmpl.cpp").render(**jinja_context)

    source_extension = "cpp" if target == Target.CPU else "cu"
    ctx.write_file(f"{class_name}.h", header)
    ctx.write_file(f"{class_name}.{source_extension}", source)
