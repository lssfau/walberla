import numpy as np
from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils import Field, FieldType
from pystencils.boundaries.boundaryhandling import create_boundary_kernel
from pystencils.boundaries.createindexlist import (
    boundary_index_array_coordinate_names, direction_member_name,
    numpy_data_type_for_boundary_object)
from pystencils.data_types import TypedSymbol, create_type
from pystencils_walberla.codegen import KernelInfo
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env


def generate_staggered_boundary(generation_context, class_name, boundary_object,
                                dim, neighbor_stencil, index_shape, target='cpu'):
    struct_name = "IndexInfo"
    boundary_object.name = class_name

    index_struct_dtype = numpy_data_type_for_boundary_object(boundary_object, dim)

    staggered_field = Field.create_generic('field', dim,
                                           np.float64 if generation_context.double_accuracy else np.float32,
                                           index_dimensions=len(index_shape), layout='c', index_shape=index_shape,
                                           field_type=FieldType.STAGGERED)

    index_field = Field('indexVector', FieldType.INDEXED, index_struct_dtype, layout=[0],
                        shape=(TypedSymbol("indexVectorSize", create_type(np.int64)), 1), strides=(1, 1))

    kernel = create_boundary_kernel(staggered_field, index_field, neighbor_stencil, boundary_object, target=target,
                                    openmp=generation_context.openmp)
    kernel.function_name = "boundary_" + boundary_object.name
    kernel.assumed_inner_stride_one = False

    # waLBerla is a 3D framework. Therefore, a zero for the z index has to be added if we work in 2D
    if dim == 2:
        stencil = ()
        for d in neighbor_stencil:
            d = d + (0,)
            stencil = stencil + (d,)
    else:
        stencil = neighbor_stencil

    stencil_info = [(i, d, ", ".join([str(e) for e in d])) for i, d in enumerate(stencil)]
    inv_dirs = []
    for direction in stencil:
        inverse_dir = tuple([-i for i in direction])
        inv_dirs.append(stencil.index(inverse_dir))

    context = {
        'class_name': boundary_object.name,
        'StructName': struct_name,
        'StructDeclaration': struct_from_numpy_dtype(struct_name, index_struct_dtype),
        'kernel': KernelInfo(kernel),
        'stencil_info': stencil_info,
        'inverse_directions': inv_dirs,
        'dim': dim,
        'target': target,
        'namespace': 'pystencils',
        'inner_or_boundary': boundary_object.inner_or_boundary
    }

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('Boundary.tmpl.h').render(**context)
    source = env.get_template('Boundary.tmpl.cpp').render(**context)

    source_extension = "cpp" if target == "cpu" else "cu"
    generation_context.write_file("{}.h".format(class_name), header)
    generation_context.write_file("{}.{}".format(class_name, source_extension), source)


def generate_staggered_flux_boundary(generation_context, class_name, boundary_object,
                                     dim, neighbor_stencil, index_shape, target='cpu'):
    struct_name = "IndexInfo"
    boundary_object.name = class_name

    index_struct_dtype = numpy_data_type_for_boundary_object(boundary_object, dim)

    staggered_field = Field.create_generic('flux', dim,
                                           np.float64 if generation_context.double_accuracy else np.float32,
                                           index_dimensions=len(index_shape), layout='c', index_shape=index_shape,
                                           field_type=FieldType.STAGGERED_FLUX)

    index_field = Field('indexVector', FieldType.INDEXED, index_struct_dtype, layout=[0],
                        shape=(TypedSymbol("indexVectorSize", create_type(np.int64)), 1), strides=(1, 1))

    kernel = create_boundary_kernel(staggered_field, index_field, neighbor_stencil, boundary_object, target=target,
                                    openmp=generation_context.openmp)
    kernel.function_name = "boundary_" + boundary_object.name
    kernel.assumed_inner_stride_one = False

    # waLBerla is a 3D framework. Therefore, a zero for the z index has to be added if we work in 2D
    if dim == 2:
        stencil = ()
        for d in neighbor_stencil:
            d = d + (0,)
            stencil = stencil + (d,)
    else:
        stencil = neighbor_stencil

    stencil_info = [(i, d, ", ".join([str(e) for e in d])) for i, d in enumerate(stencil)]
    inv_dirs = []
    for direction in stencil:
        inverse_dir = tuple([-i for i in direction])
        inv_dirs.append(stencil.index(inverse_dir))

    context = {
        'class_name': boundary_object.name,
        'StructName': struct_name,
        'StructDeclaration': struct_from_numpy_dtype(struct_name, index_struct_dtype),
        'kernel': KernelInfo(kernel),
        'stencil_info': stencil_info,
        'inverse_directions': inv_dirs,
        'dim': dim,
        'target': target,
        'namespace': 'pystencils',
        'inner_or_boundary': boundary_object.inner_or_boundary
    }

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('Boundary.tmpl.h').render(**context)
    source = env.get_template('Boundary.tmpl.cpp').render(**context)

    source_extension = "cpp" if target == "cpu" else "cu"
    generation_context.write_file("{}.h".format(class_name), header)
    generation_context.write_file("{}.{}".format(class_name, source_extension), source)


def struct_from_numpy_dtype(struct_name, numpy_dtype):
    result = "struct %s { \n" % (struct_name,)

    equality_compare = []
    constructor_params = []
    constructor_initializer_list = []
    for name, (sub_type, offset) in numpy_dtype.fields.items():
        pystencils_type = create_type(sub_type)
        result += "    %s %s;\n" % (pystencils_type, name)
        if name in boundary_index_array_coordinate_names or name == direction_member_name:
            constructor_params.append("%s %s_" % (pystencils_type, name))
            constructor_initializer_list.append("%s(%s_)" % (name, name))
        else:
            constructor_initializer_list.append("%s()" % name)
        if pystencils_type.is_float():
            equality_compare.append("floatIsEqual(%s, o.%s)" % (name, name))
        else:
            equality_compare.append("%s == o.%s" % (name, name))

    result += "    %s(%s) : %s {}\n" % \
              (struct_name, ", ".join(constructor_params), ", ".join(constructor_initializer_list))
    result += "    bool operator==(const %s & o) const {\n        return %s;\n    }\n" % \
              (struct_name, " && ".join(equality_compare))
    result += "};\n"
    return result
