from typing import Sequence, Union


from pystencils import Target, Assignment, AssignmentCollection
from pystencils import create_kernel, create_staggered_kernel

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.kernel_selection import KernelCallNode, KernelFamily, HighLevelInterfaceSpec
from pystencils_walberla.utility import config_from_context


def function_generator(ctx: CodeGenerationContext, class_name: str,
                       assignments: Union[Sequence[Assignment], AssignmentCollection],
                       namespace: str = 'pystencils', staggered=False, field_swaps=None, varying_parameters=(),
                       ghost_layers_to_include=0,
                       target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None,
                       max_threads=None,
                       **create_kernel_params):
    return lambda: __function_generator(ctx, class_name, assignments,
                                        namespace, staggered, field_swaps, varying_parameters,
                                        ghost_layers_to_include,
                                        target, data_type, cpu_openmp, cpu_vectorize_info, max_threads,
                                        **create_kernel_params)


def __function_generator(ctx: CodeGenerationContext, class_name: str,
                         assignments: Union[Sequence[Assignment], AssignmentCollection],
                         namespace: str = 'pystencils', staggered=False, field_swaps=None, varying_parameters=(),
                         ghost_layers_to_include=0,
                         target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None,
                         max_threads=None,
                         **create_kernel_params):
    if staggered:
        assert 'omp_single_loop' not in create_kernel_params

    create_kernel_params['omp_single_loop'] = False
    config = config_from_context(ctx, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 cpu_vectorize_info=cpu_vectorize_info, **create_kernel_params)

    if not staggered:
        ast = create_kernel(assignments, config=config)
    else:
        # This should not be necessary but create_staggered_kernel does not take a config at the moment ...
        ast = create_staggered_kernel(assignments, **config.__dict__)

    ast.function_name = class_name.lower()

    all_field_names = [f.name for f in ast.fields_accessed]
    all_field_names.sort()

    temporary_fields = [f for f in all_field_names if "_tmp" in f]

    if field_swaps is None:
        field_swaps = []
        for field_name in all_field_names:
            if field_name + "_tmp" in temporary_fields:
                field_swaps.append((field_name, field_name + "_tmp"))

    selection_tree = KernelCallNode(ast)
    kernel_family = KernelFamily(selection_tree, class_name,
                                 temporary_fields, field_swaps, varying_parameters)

    representative_field = {p.field_name for p in kernel_family.parameters if p.is_field_parameter}
    representative_field = sorted(representative_field)[0]

    interface_spec = HighLevelInterfaceSpec(kernel_family.kernel_selection_parameters, ())

    jinja_context = {
        'kernel': kernel_family,
        'namespace': namespace,
        'function_name': class_name,
        'field': representative_field,
        'ghost_layers_to_include': ghost_layers_to_include,
        'interface_spec': interface_spec,
        'max_threads': max_threads
    }
    return jinja_context
