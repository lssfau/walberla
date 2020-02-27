import jinja2
import sympy as sp

from pystencils import TypedSymbol
from pystencils.backends.cbackend import generate_c
from pystencils.backends.cuda_backend import CudaSympyPrinter
from pystencils.data_types import get_base_type
from pystencils.field import FieldType
from pystencils.kernelparameters import SHAPE_DTYPE
from pystencils.sympyextensions import prod

temporary_fieldMemberTemplate = """
private: std::set< {type} *, field::SwapableCompare< {type} * > > cache_{original_field_name}_;"""

temporary_fieldTemplate = """
// Getting temporary field {tmp_field_name}
auto it = cache_{original_field_name}_.find( {original_field_name} );
if( it != cache_{original_field_name}_.end() )
{{
    {tmp_field_name} = *it;
}}
else 
{{
    {tmp_field_name} = {original_field_name}->cloneUninitialized();
    cache_{original_field_name}_.insert({tmp_field_name});
}}
"""

temporary_constructor = """
~{class_name}() {{  {contents} }}
"""

delete_loop = """
    for(auto p: cache_{original_field_name}_) {{
        delete p;
    }}
"""


def make_field_type(dtype, f_size, is_gpu):
    if is_gpu:
        return "cuda::GPUField<%s>" % (dtype,)
    else:
        return "GhostLayerField<%s, %d>" % (dtype, f_size)


def get_field_fsize(field):
    """Determines the size of the index coordinate. Since walberla fields only support one index dimension,
    pystencils fields with multiple index dimensions are linearized to a single index dimension.
    """
    assert field.has_fixed_index_shape, \
        "All Fields have to be created with fixed index coordinate shape using index_shape=(q,) " + str(field.name)

    if field.index_dimensions == 0:
        return 1
    else:
        return prod(field.index_shape)


def get_field_stride(param):
    field = param.fields[0]
    type_str = get_base_type(param.symbol.dtype).base_name
    stride_names = ['xStride()', 'yStride()', 'zStride()', 'fStride()']
    stride_names = ["%s(%s->%s)" % (type_str, param.field_name, e) for e in stride_names]
    strides = stride_names[:field.spatial_dimensions]
    if field.index_dimensions > 0:
        additional_strides = [1]
        for shape in reversed(field.index_shape[1:]):
            additional_strides.append(additional_strides[-1] * shape)
        assert len(additional_strides) == field.index_dimensions
        f_stride_name = stride_names[-1]
        strides.extend(["%s(%d * %s)" % (type_str, e, f_stride_name) for e in reversed(additional_strides)])
    return strides[param.symbol.coordinate]


def generate_declaration(kernel_info, target='cpu'):
    """Generates the declaration of the kernel function"""
    ast = kernel_info.ast
    result = generate_c(ast, signature_only=True, dialect='cuda' if target == 'gpu' else 'c') + ";"
    result = "namespace internal_%s {\n%s\n}" % (ast.function_name, result,)
    return result


def generate_definition(kernel_info, target='cpu'):
    """Generates the definition (i.e. implementation) of the kernel function"""
    ast = kernel_info.ast
    result = generate_c(ast, dialect='cuda' if target == 'gpu' else 'c')
    result = "namespace internal_%s {\nstatic %s\n}" % (ast.function_name, result)
    return result


def field_extraction_code(field, is_temporary, declaration_only=False,
                          no_declaration=False, is_gpu=False):
    """Returns code string for getting a field pointer.

    This can happen in two ways: either the field is extracted from a walberla block, or a temporary field to swap is
    created.

    Args:
        field: the field for which the code should be created
        is_temporary: new_filtered field from block (False) or create a temporary copy of an existing field (True)
        declaration_only: only create declaration instead of the full code
        no_declaration: create the extraction code, and assume that declarations are elsewhere
        is_gpu: if the field is a GhostLayerField or a GpuField
    """

    # Determine size of f coordinate which is a template parameter
    f_size = get_field_fsize(field)
    field_name = field.name
    dtype = get_base_type(field.dtype)
    field_type = make_field_type(dtype, f_size, is_gpu)

    if not is_temporary:
        dtype = get_base_type(field.dtype)
        field_type = make_field_type(dtype, f_size, is_gpu)
        if declaration_only:
            return "%s * %s;" % (field_type, field_name)
        else:
            prefix = "" if no_declaration else "auto "
            return "%s%s = block->uncheckedFastGetData< %s >(%sID);" % (prefix, field_name, field_type, field_name)
    else:
        assert field_name.endswith('_tmp')
        original_field_name = field_name[:-len('_tmp')]
        if declaration_only:
            return "%s * %s;" % (field_type, field_name)
        else:
            declaration = "{type} * {tmp_field_name};".format(type=field_type, tmp_field_name=field_name)
            tmp_field_str = temporary_fieldTemplate.format(original_field_name=original_field_name,
                                                           tmp_field_name=field_name, type=field_type)
            return tmp_field_str if no_declaration else declaration + tmp_field_str


@jinja2.contextfilter
def generate_block_data_to_field_extraction(ctx, kernel_info, parameters_to_ignore=(), parameters=None,
                                            declarations_only=False, no_declarations=False):
    """Generates code that extracts all required fields of a kernel from a walberla block storage."""
    if parameters is not None:
        assert parameters_to_ignore == ()
        field_parameters = []
        for param in kernel_info.parameters:
            if param.is_field_pointer and param.field_name in parameters:
                field_parameters.append(param.fields[0])
    else:
        field_parameters = []
        for param in kernel_info.parameters:
            if param.is_field_pointer and param.field_name not in parameters_to_ignore:
                field_parameters.append(param.fields[0])

    normal_fields = {f for f in field_parameters if f.name not in kernel_info.temporary_fields}
    temporary_fields = {f for f in field_parameters if f.name in kernel_info.temporary_fields}

    args = {
        'declaration_only': declarations_only,
        'no_declaration': no_declarations,
        'is_gpu': ctx['target'] == 'gpu',
    }
    result = "\n".join(field_extraction_code(field=field, is_temporary=False, **args) for field in normal_fields) + "\n"
    result += "\n".join(field_extraction_code(field=field, is_temporary=True, **args) for field in temporary_fields)
    return result


def generate_refs_for_kernel_parameters(kernel_info, prefix, parameters_to_ignore):
    symbols = {p.field_name for p in kernel_info.parameters if p.is_field_pointer}
    symbols.update(p.symbol.name for p in kernel_info.parameters if not p.is_field_parameter)
    symbols.difference_update(parameters_to_ignore)
    return "\n".join("auto & %s = %s%s;" % (s, prefix, s) for s in symbols)


@jinja2.contextfilter
def generate_call(ctx, kernel_info, ghost_layers_to_include=0, cell_interval=None, stream='0',
                  spatial_shape_symbols=()):
    """Generates the function call to a pystencils kernel

    Args:
        kernel_info:
        ghost_layers_to_include: if left to 0, only the inner part of the ghost layer field is looped over
                                 a CHECK is inserted that the field has as many ghost layers as the pystencils AST
                                 needs. This parameter specifies how many ghost layers the kernel should view as
                                 "inner area". The ghost layer field has to have the required number of ghost layers
                                 remaining. Parameter has to be left to default if cell_interval is given.
        cell_interval: Defines the name (string) of a walberla CellInterval object in scope,
                       that defines the inner region for the kernel to loop over. Parameter has to be left to default
                       if ghost_layers_to_include is specified.
        stream: optional name of cuda stream variable
        spatial_shape_symbols: relevant only for gpu kernels - to determine CUDA block and grid sizes the iteration
                               region (i.e. field shape) has to be known. This can normally be inferred by the kernel
                               parameters - however in special cases like boundary conditions a manual specification
                               may be necessary.
    """
    assert isinstance(ghost_layers_to_include, str) or ghost_layers_to_include >= 0
    ast = kernel_info.ast
    ast_params = kernel_info.parameters
    is_cpu = ctx['target'] == 'cpu'

    ghost_layers_to_include = sp.sympify(ghost_layers_to_include)
    if ast.ghost_layers is None:
        required_ghost_layers = 0
    else:
        # ghost layer info is ((x_gl_front, x_gl_end), (y_gl_front, y_gl_end).. )
        required_ghost_layers = max(max(ast.ghost_layers))

    kernel_call_lines = []

    def get_start_coordinates(field_object):
        if cell_interval is None:
            return [-ghost_layers_to_include - required_ghost_layers] * field_object.spatial_dimensions
        else:
            assert ghost_layers_to_include == 0
            return [sp.Symbol("{ci}.{coord}Min()".format(coord=coord_name, ci=cell_interval)) - required_ghost_layers
                    for coord_name in ('x', 'y', 'z')]

    def get_end_coordinates(field_object):
        if cell_interval is None:
            shape_names = ['xSize()', 'ySize()', 'zSize()'][:field_object.spatial_dimensions]
            offset = 2 * ghost_layers_to_include + 2 * required_ghost_layers
            return ["cell_idx_c(%s->%s) + %s" % (field_object.name, e, offset) for e in shape_names]
        else:
            assert ghost_layers_to_include == 0
            return ["cell_idx_c({ci}.{coord}Size()) + {gl}".format(coord=coord_name, ci=cell_interval,
                                                                   gl=2 * required_ghost_layers)
                    for coord_name in ('x', 'y', 'z')]

    for param in ast_params:
        if param.is_field_parameter and FieldType.is_indexed(param.fields[0]):
            continue

        if param.is_field_pointer:
            field = param.fields[0]
            if field.field_type == FieldType.BUFFER:
                kernel_call_lines.append("%s %s = %s;" % (param.symbol.dtype, param.symbol.name, param.field_name))
            else:
                coordinates = get_start_coordinates(field)
                actual_gls = "int_c(%s->nrOfGhostLayers())" % (param.field_name, )
                coord_set = set(coordinates)
                coord_set = sorted(coord_set, key=lambda e: str(e))
                for c in coord_set:
                    kernel_call_lines.append("WALBERLA_ASSERT_GREATER_EQUAL(%s, -%s);" %
                                             (c, actual_gls))
                while len(coordinates) < 4:
                    coordinates.append(0)
                coordinates = tuple(coordinates)
                kernel_call_lines.append("%s %s = %s->dataAt(%s, %s, %s, %s);" %
                                         ((param.symbol.dtype, param.symbol.name, param.field_name) + coordinates))
        elif param.is_field_stride:
            casted_stride = get_field_stride(param)
            type_str = param.symbol.dtype.base_name
            kernel_call_lines.append("const %s %s = %s;" % (type_str, param.symbol.name, casted_stride))
        elif param.is_field_shape:
            coord = param.symbol.coordinate
            field = param.fields[0]
            type_str = param.symbol.dtype.base_name
            shape = "%s(%s)" % (type_str, get_end_coordinates(field)[coord])
            assert coord < 3
            max_value = "%s->%sSizeWithGhostLayer()" % (field.name, ('x', 'y', 'z')[coord])
            kernel_call_lines.append("WALBERLA_ASSERT_GREATER_EQUAL(%s, %s);" % (max_value, shape))
            kernel_call_lines.append("const %s %s = %s;" % (type_str, param.symbol.name, shape))

    call_parameters = ", ".join([p.symbol.name for p in ast_params])
    if not is_cpu:
        if not spatial_shape_symbols:
            spatial_shape_symbols = [p.symbol for p in ast_params if p.is_field_shape]
            spatial_shape_symbols.sort(key=lambda e: e.coordinate)
        else:
            spatial_shape_symbols = [TypedSymbol(s, SHAPE_DTYPE) for s in spatial_shape_symbols]

        indexing_dict = ast.indexing.call_parameters(spatial_shape_symbols)
        sp_printer_c = CudaSympyPrinter()
        kernel_call_lines += [
            "dim3 _block(int(%s), int(%s), int(%s));" % tuple(sp_printer_c.doprint(e) for e in indexing_dict['block']),
            "dim3 _grid(int(%s), int(%s), int(%s));" % tuple(sp_printer_c.doprint(e) for e in indexing_dict['grid']),
            "internal_%s::%s<<<_grid, _block, 0, %s>>>(%s);" % (ast.function_name, ast.function_name,
                                                                stream, call_parameters),
        ]
    else:
        kernel_call_lines.append("internal_%s::%s(%s);" % (ast.function_name, ast.function_name, call_parameters))
    return "\n".join(kernel_call_lines)


def generate_swaps(kernel_info):
    """Generates code to swap main fields with temporary fields"""
    swaps = ""
    for src, dst in kernel_info.field_swaps:
        swaps += "%s->swapDataPointers(%s);\n" % (src, dst)
    return swaps


def generate_constructor_initializer_list(kernel_info, parameters_to_ignore=None):
    if parameters_to_ignore is None:
        parameters_to_ignore = []

    parameters_to_ignore += kernel_info.temporary_fields

    parameter_initializer_list = []
    for param in kernel_info.parameters:
        if param.is_field_pointer and param.field_name not in parameters_to_ignore:
            parameter_initializer_list.append("%sID(%sID_)" % (param.field_name, param.field_name))
        elif not param.is_field_parameter and param.symbol.name not in parameters_to_ignore:
            parameter_initializer_list.append("%s(%s_)" % (param.symbol.name, param.symbol.name))
    return ", ".join(parameter_initializer_list)


def generate_constructor_parameters(kernel_info, parameters_to_ignore=None):
    if parameters_to_ignore is None:
        parameters_to_ignore = []
        
    varying_parameters = []
    if hasattr(kernel_info, 'varying_parameters'):
        varying_parameters = kernel_info.varying_parameters
    varying_parameter_names = tuple(e[1] for e in varying_parameters)
    parameters_to_ignore += kernel_info.temporary_fields + varying_parameter_names

    parameter_list = []
    for param in kernel_info.parameters:
        if param.is_field_pointer and param.field_name not in parameters_to_ignore:
            parameter_list.append("BlockDataID %sID_" % (param.field_name, ))
        elif not param.is_field_parameter and param.symbol.name not in parameters_to_ignore:
            parameter_list.append("%s %s_" % (param.symbol.dtype, param.symbol.name,))
    varying_parameters = ["%s %s_" % e for e in varying_parameters]
    return ", ".join(parameter_list + varying_parameters)


@jinja2.contextfilter
def generate_members(ctx, kernel_info, parameters_to_ignore=(), only_fields=False):
    ast = kernel_info.ast
    fields = {f.name: f for f in ast.fields_accessed}

    params_to_skip = tuple(parameters_to_ignore) + tuple(kernel_info.temporary_fields)
    params_to_skip += tuple(e[1] for e in kernel_info.varying_parameters)
    is_gpu = ctx['target'] == 'gpu'

    result = []
    for param in kernel_info.parameters:
        if only_fields and not param.is_field_parameter:
            continue
        if param.is_field_pointer and param.field_name not in params_to_skip:
            result.append("BlockDataID %sID;" % (param.field_name, ))
        elif not param.is_field_parameter and param.symbol.name not in params_to_skip:
            result.append("%s %s;" % (param.symbol.dtype, param.symbol.name,))

    for field_name in kernel_info.temporary_fields:
        f = fields[field_name]
        if field_name in parameters_to_ignore:
            continue
        assert field_name.endswith('_tmp')
        original_field_name = field_name[:-len('_tmp')]
        f_size = get_field_fsize(f)
        field_type = make_field_type(get_base_type(f.dtype), f_size, is_gpu)
        result.append(temporary_fieldMemberTemplate.format(type=field_type, original_field_name=original_field_name))

    if hasattr(kernel_info, 'varying_parameters'):
        result.extend(["%s %s;" % e for e in kernel_info.varying_parameters])

    return "\n".join(result)


def generate_destructor(kernel_info, class_name):
    if not kernel_info.temporary_fields:
        return ""
    else:
        contents = ""
        for field_name in kernel_info.temporary_fields:
            contents += delete_loop.format(original_field_name=field_name[:-len('_tmp')])
        return temporary_constructor.format(contents=contents, class_name=class_name)


def add_pystencils_filters_to_jinja_env(jinja_env):
    jinja_env.filters['generate_definition'] = generate_definition
    jinja_env.filters['generate_declaration'] = generate_declaration
    jinja_env.filters['generate_members'] = generate_members
    jinja_env.filters['generate_constructor_parameters'] = generate_constructor_parameters
    jinja_env.filters['generate_constructor_initializer_list'] = generate_constructor_initializer_list
    jinja_env.filters['generate_call'] = generate_call
    jinja_env.filters['generate_block_data_to_field_extraction'] = generate_block_data_to_field_extraction
    jinja_env.filters['generate_swaps'] = generate_swaps
    jinja_env.filters['generate_refs_for_kernel_parameters'] = generate_refs_for_kernel_parameters
    jinja_env.filters['generate_destructor'] = generate_destructor
