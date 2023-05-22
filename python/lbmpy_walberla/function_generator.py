from pystencils_walberla.kernel_selection import KernelCallNode, KernelFamily, HighLevelInterfaceSpec


def kernel_family_function_generator(class_name: str, kernel_family: KernelFamily,
                                     namespace: str = 'lbm', max_threads: int = None):

    return lambda: __function_generator(class_name, kernel_family, namespace, max_threads)


def __function_generator(class_name: str, kernel_family: KernelFamily,
                         namespace: str = 'lbm', max_threads: int = None):

    representative_field = {p.field_name for p in kernel_family.parameters if p.is_field_parameter}
    representative_field = sorted(representative_field)[0]

    interface_spec = HighLevelInterfaceSpec(kernel_family.kernel_selection_parameters, ())

    jinja_context = {
        'kernel': kernel_family,
        'namespace': namespace,
        'function_name': class_name,
        'field': representative_field,
        'interface_spec': interface_spec,
        'max_threads': max_threads
    }
    return jinja_context
