import numpy as np
from jinja2 import Environment, PackageLoader, StrictUndefined
from pystencils import Field, FieldType, Target
from pystencils.boundaries.boundaryhandling import create_boundary_kernel
from pystencils.boundaries.createindexlist import numpy_data_type_for_boundary_object
from pystencils.typing import TypedSymbol, create_type

from pystencils_walberla.utility import config_from_context, struct_from_numpy_dtype
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.additional_data_handler import AdditionalDataHandler
from pystencils_walberla.kernel_selection import (
    KernelFamily, AbstractKernelSelectionNode, KernelCallNode, HighLevelInterfaceSpec)
from pystencils.astnodes import KernelFunction


def generate_boundary(generation_context,
                      class_name,
                      boundary_object,
                      field_name,
                      neighbor_stencil,
                      index_shape,
                      field_type=FieldType.GENERIC,
                      kernel_creation_function=None,
                      target=Target.CPU,
                      data_type=None,
                      cpu_openmp=None,
                      namespace='pystencils',
                      additional_data_handler=None,
                      interface_mappings=(),
                      generate_functor=True,
                      layout='fzyx',
                      field_timestep=None,
                      **create_kernel_params):

    if boundary_object.additional_data and additional_data_handler is None:
        raise ValueError("Boundary object has additional data but you have not provided an AdditionalDataHandler.")

    struct_name = "IndexInfo"
    boundary_object.name = class_name
    dim = neighbor_stencil.D

    config = config_from_context(generation_context, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 **create_kernel_params)
    create_kernel_params = config.__dict__
    del create_kernel_params['target']
    del create_kernel_params['index_fields']
    del create_kernel_params['default_number_int']
    del create_kernel_params['skip_independence_check']

    field_data_type = config.data_type[field_name].numpy_dtype

    index_struct_dtype = numpy_data_type_for_boundary_object(boundary_object, dim)

    field = Field.create_generic(field_name, dim,
                                 field_data_type,
                                 index_dimensions=len(index_shape), layout=layout, index_shape=index_shape,
                                 field_type=field_type)

    index_field = Field('indexVector', FieldType.INDEXED, index_struct_dtype, layout=[0],
                        shape=(TypedSymbol("indexVectorSize", create_type("int32")), 1), strides=(1, 1))

    if not kernel_creation_function:
        kernel_creation_function = create_boundary_kernel

    kernel = kernel_creation_function(field, index_field, neighbor_stencil, boundary_object,
                                      target=target, **create_kernel_params)

    if isinstance(kernel, KernelFunction):
        kernel.function_name = f"boundary_{boundary_object.name}"
        selection_tree = KernelCallNode(kernel)
    elif isinstance(kernel, AbstractKernelSelectionNode):
        selection_tree = kernel
    else:
        raise ValueError(f"kernel_creation_function returned wrong type: {kernel.__class__}")

    kernel_family = KernelFamily(selection_tree, class_name, field_timestep=field_timestep)
    selection_parameters = kernel_family.kernel_selection_parameters if field_timestep is None else []
    interface_spec = HighLevelInterfaceSpec(selection_parameters, interface_mappings)

    if additional_data_handler is None:
        additional_data_handler = AdditionalDataHandler(stencil=neighbor_stencil)

    default_dtype = config.data_type.default_factory()
    is_float = True if issubclass(default_dtype.numpy_dtype.type, np.float32) else False

    context = {
        'kernel': kernel_family,
        'class_name': boundary_object.name,
        'interface_spec': interface_spec,
        'generate_functor': generate_functor,
        'StructName': struct_name,
        'StructDeclaration': struct_from_numpy_dtype(struct_name, index_struct_dtype),
        'dim': dim,
        'target': target.name.lower(),
        'namespace': namespace,
        'inner_or_boundary': boundary_object.inner_or_boundary,
        'single_link': boundary_object.single_link,
        'additional_data_handler': additional_data_handler,
        'dtype': "double" if is_float else "float",
        'layout': layout,
        'index_shape': index_shape
    }

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('Boundary.tmpl.h').render(**context)
    source = env.get_template('Boundary.tmpl.cpp').render(**context)

    source_extension = "cu" if target == Target.GPU and generation_context.cuda else "cpp"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)

    return context


def generate_staggered_boundary(generation_context, class_name, boundary_object,
                                dim, neighbor_stencil, index_shape, target=Target.CPU, **kwargs):
    assert dim == len(neighbor_stencil[0])
    generate_boundary(generation_context, class_name, boundary_object, 'field', neighbor_stencil, index_shape,
                      FieldType.STAGGERED, target=target, **kwargs)


def generate_staggered_flux_boundary(generation_context, class_name, boundary_object,
                                     dim, neighbor_stencil, index_shape, target=Target.CPU, **kwargs):
    assert dim == len(neighbor_stencil[0])
    generate_boundary(generation_context, class_name, boundary_object, 'flux', neighbor_stencil, index_shape,
                      FieldType.STAGGERED_FLUX, target=target, **kwargs)


