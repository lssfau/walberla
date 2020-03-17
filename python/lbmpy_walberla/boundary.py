import numpy as np
from jinja2 import Environment, PackageLoader, StrictUndefined

from lbmpy.boundaries.boundaryhandling import create_lattice_boltzmann_boundary_kernel
from lbmpy_walberla.walberla_lbm_generation import KernelInfo
from pystencils import Field, FieldType
from pystencils.boundaries.createindexlist import (
    boundary_index_array_coordinate_names, direction_member_name,
    numpy_data_type_for_boundary_object)
from pystencils.data_types import TypedSymbol, create_type
from pystencils_walberla.codegen import default_create_kernel_parameters
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env


def generate_boundary(generation_context, class_name, boundary_object, lb_method, **create_kernel_params):
    struct_name = "IndexInfo"
    boundary_object.name = class_name

    create_kernel_params = default_create_kernel_parameters(generation_context, create_kernel_params)
    target = create_kernel_params['target']

    index_struct_dtype = numpy_data_type_for_boundary_object(boundary_object, lb_method.dim)

    pdf_field = Field.create_generic('pdfs', lb_method.dim,
                                     np.float64 if generation_context.double_accuracy else np.float32,
                                     index_dimensions=1, layout='fzyx', index_shape=[len(lb_method.stencil)])

    index_field = Field('indexVector', FieldType.INDEXED, index_struct_dtype, layout=[0],
                        shape=(TypedSymbol("indexVectorSize", create_type(np.int64)), 1), strides=(1, 1))

    kernel = create_lattice_boltzmann_boundary_kernel(pdf_field, index_field, lb_method, boundary_object, target=target,
                                                      openmp=generation_context.openmp)
    kernel.function_name = "boundary_" + boundary_object.name

    # waLBerla is a 3D framework. Therefore, a zero for the z index has to be added if we work in 2D
    if lb_method.dim == 2:
        stencil = ()
        for d in lb_method.stencil:
            d = d + (0,)
            stencil = stencil + (d,)
    else:
        stencil = lb_method.stencil

    stencil_info = [(i, ", ".join([str(e) for e in d])) for i, d in enumerate(stencil)]

    context = {
        'class_name': boundary_object.name,
        'StructName': struct_name,
        'StructDeclaration': struct_from_numpy_dtype(struct_name, index_struct_dtype),
        'kernel': KernelInfo(kernel),
        'stencil_info': stencil_info,
        'dim': lb_method.dim,
        'target': target,
        'namespace': 'lbm',
    }

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('Boundary.tmpl.h').render(**context)
    source = env.get_template('Boundary.tmpl.cpp').render(**context)

    source_extension = "cpp" if create_kernel_params.get("target", "cpu") == "cpu" else "cu"
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
