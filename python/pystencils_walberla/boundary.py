import numpy as np
from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils_walberla.codegen import KernelInfo
from pystencils import Field, FieldType
from pystencils.boundaries.createindexlist import (
    boundary_index_array_coordinate_names, direction_member_name,
    numpy_data_type_for_boundary_object)
from pystencils.data_types import TypedSymbol, create_type
from pystencils_walberla.codegen import default_create_kernel_parameters
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env


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
