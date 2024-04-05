# For backward compatibility with version < 3.0.0
try:
    from jinja2 import pass_context as jinja2_context_decorator
except ImportError:
    from jinja2 import contextfilter as jinja2_context_decorator

from collections.abc import Iterable
import sympy as sp

from pystencils import Target, Backend
from pystencils.backends.cbackend import generate_c
from pystencils.typing import TypedSymbol, get_base_type
from pystencils.field import FieldType
from pystencils.sympyextensions import prod

temporary_fieldPointerTemplate = """{type}"""

temporary_fieldMemberTemplate = """
private: std::set< {type} *, field::SwapableCompare< {type} * > > cache_{original_field_name}_;"""

temporary_fieldTemplate = """
{{
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

standard_parameter_registration = """
for (uint_t level = 0; level < blocks->getNumberOfLevels(); level++)
{{
    const {dtype} level_scale_factor = {dtype}(uint_t(1) << level);
    const {dtype} one                = {dtype}(1.0);
    const {dtype} half               = {dtype}(0.5);
    
    {name}Vector.push_back( {dtype}({name} / (level_scale_factor * (-{name} * half + one) + {name} * half)) );
}}
"""


# the target will enter the jinja filters as string. The reason for that is, that is not easy to work with the
# enum in the template files.
def translate_target(target):
    if isinstance(target, Target):
        return target
    elif isinstance(target, str):
        return Target[target.upper()]
    else:
        raise ValueError(f"The type of the target {type(target)} is not supported")


def make_field_type(dtype, f_size, is_gpu):
    if is_gpu:
        return f"gpu::GPUField<{dtype}>"
    else:
        return f"field::GhostLayerField<{dtype}, {f_size}>"


def field_type(field, is_gpu=False):
    dtype = get_base_type(field.dtype)
    f_size = get_field_fsize(field)
    return make_field_type(dtype, f_size, is_gpu)


def get_field_fsize(field):
    """Determines the size of the index coordinate. Since walberla fields only support one index dimension,
    pystencils fields with multiple index dimensions are linearized to a single index dimension.
    """
    assert field.has_fixed_index_shape, \
        f"All Fields have to be created with fixed index coordinate shape using index_shape=(q,) {str(field.name)}"

    if field.index_dimensions == 0:
        return 1
    else:
        return prod(field.index_shape)


def get_field_stride(param):
    field = param.fields[0]
    type_str = get_base_type(param.symbol.dtype).c_name
    stride_names = ['xStride()', 'yStride()', 'zStride()', 'fStride()']
    stride_names = [f"{type_str}({param.field_name}->{e})" for e in stride_names]
    strides = stride_names[:field.spatial_dimensions]
    if field.index_dimensions > 0:
        additional_strides = [1]
        for shape in reversed(field.index_shape[1:]):
            additional_strides.append(additional_strides[-1] * shape)
        assert len(additional_strides) == field.index_dimensions
        f_stride_name = stride_names[-1]
        strides.extend([f"{type_str}({e} * {f_stride_name})" for e in reversed(additional_strides)])
    return strides[param.symbol.coordinate]


def generate_declaration(kernel_info, target=Target.CPU):
    """Generates the declaration of the kernel function"""
    target = translate_target(target)
    ast = kernel_info.ast
    result = generate_c(ast, signature_only=True, dialect=Backend.CUDA if target == Target.GPU else Backend.C) + ";"
    result = "namespace internal_%s {\n%s\n}" % (ast.function_name, result,)
    return result


def generate_definition(kernel_info, target=Target.CPU):
    """Generates the definition (i.e. implementation) of the kernel function"""
    target = translate_target(target)
    ast = kernel_info.ast
    result = generate_c(ast, dialect=Backend.CUDA if target == Target.GPU else Backend.C)
    result = "namespace internal_%s {\nstatic %s\n}" % (ast.function_name, result)
    return result


def generate_declarations(kernel_family, target=Target.CPU):
    target = translate_target(target)
    declarations = []
    for ast in kernel_family.all_asts:
        code = generate_c(ast, signature_only=True, dialect=Backend.CUDA if target == Target.GPU else Backend.C) + ";"
        code = "namespace internal_%s {\n%s\n}\n" % (ast.function_name, code,)
        declarations.append(code)
    return "\n".join(declarations)


def generate_definitions(kernel_family, target=Target.CPU, max_threads=None):
    target = translate_target(target)
    definitions = []
    for ast in kernel_family.all_asts:
        code = generate_c(ast, dialect=Backend.CUDA if target == Target.GPU else Backend.C)
        if max_threads is not None and target == Target.GPU:
            assert isinstance(max_threads, int), "maximal number of threads should be an integer value"
            index = code.find('FUNC_PREFIX') + len("FUNC_PREFIX ")
            code = code[:index] + f'__launch_bounds__({max_threads}) ' + code[index:]
        code = "namespace internal_%s {\nstatic %s\n}\n" % (ast.function_name, code)
        definitions.append(code)
    return "\n".join(definitions)


def field_extraction_code(field, is_temporary, declaration_only=False,
                          no_declaration=False, is_gpu=False, update_member=False):
    """Returns code string for getting a field pointer.

    This can happen in two ways: either the field is extracted from a walberla block, or a temporary field to swap is
    created.

    Args:
        field: the field for which the code should be created
        is_temporary: new_filtered field from block (False) or create a temporary copy of an existing field (True)
        declaration_only: only create declaration instead of the full code
        no_declaration: create the extraction code, and assume that declarations are elsewhere
        is_gpu: if the field is a GhostLayerField or a GpuField
        update_member: specify if function is used inside a constructor; add _ to members
    """
    wlb_field_type = field_type(field, is_gpu)

    if not is_temporary:
        if declaration_only:
            return f"{wlb_field_type} * {field.name}_;"
        else:
            prefix = "" if no_declaration else "auto "
            if update_member:
                return f"{prefix}{field.name}_ = block->getData< {wlb_field_type} >({field.name}ID);"
            else:
                return f"{prefix}{field.name} = block->getData< {wlb_field_type} >({field.name}ID);"
    else:
        assert field.name.endswith('_tmp')
        original_field_name = field.name[:-len('_tmp')]
        if declaration_only:
            return f"{wlb_field_type} * {field.name}_;"
        else:
            declaration = f"{wlb_field_type} * {field.name};"
            tmp_field_str = temporary_fieldTemplate.format(original_field_name=original_field_name,
                                                           tmp_field_name=field.name, type=wlb_field_type)
            return tmp_field_str if no_declaration else declaration + tmp_field_str


# TODO fields are not sorted
@jinja2_context_decorator
def generate_block_data_to_field_extraction(ctx, kernel_info, parameters_to_ignore=(), parameters=None,
                                            declarations_only=False, no_declarations=False, update_member=False):
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

    target = translate_target(ctx['target'])

    args = {
        'declaration_only': declarations_only,
        'no_declaration': no_declarations,
        'is_gpu': target == Target.GPU,
    }
    result = "\n".join(
        field_extraction_code(field=field, is_temporary=False, update_member=update_member, **args) for field in
        normal_fields) + "\n"
    result += "\n".join(
        field_extraction_code(field=field, is_temporary=True, update_member=update_member, **args) for field in
        temporary_fields)
    return result


def generate_refs_for_kernel_parameters(kernel_info, prefix, parameters_to_ignore=(), ignore_fields=False,
                                        parameter_registration=None):
    symbols = {p.field_name for p in kernel_info.parameters if p.is_field_pointer and not ignore_fields}
    symbols.update(p.symbol.name for p in kernel_info.parameters if not p.is_field_parameter)
    symbols.difference_update(parameters_to_ignore)
    type_information = {p.symbol.name: p.symbol.dtype for p in kernel_info.parameters if not p.is_field_parameter}
    result = []
    registered_parameters = [] if not parameter_registration else parameter_registration.scaling_info
    for s in symbols:
        if s in registered_parameters:
            dtype = type_information[s].c_name
            result.append("const uint_t level = block->getBlockStorage().getLevel(*block);")
            result.append(f"{dtype} & {s} = {s}Vector[level];")
        else:
            result.append(f"auto & {s} = {prefix}{s}_;")
    return "\n".join(result)


@jinja2_context_decorator
def generate_call(ctx, kernel, ghost_layers_to_include=0, cell_interval=None, stream='0',
                  spatial_shape_symbols=()):
    """Generates the function call to a pystencils kernel

    Args:
        ctx: code generation context
        kernel: pystencils kernel
        ghost_layers_to_include: if left to 0, only the inner part of the ghost layer field is looped over
                                 a CHECK is inserted that the field has as many ghost layers as the pystencils AST
                                 needs. This parameter specifies how many ghost layers the kernel should view as
                                 "inner area". The ghost layer field has to have the required number of ghost layers
                                 remaining. Parameter has to be left to default if cell_interval is given.
        cell_interval: Defines the name (string) of a walberla CellInterval object in scope,
                       that defines the inner region for the kernel to loop over. Parameter has to be left to default
                       if ghost_layers_to_include is specified.
        stream: optional name of gpu stream variable
        spatial_shape_symbols: relevant only for gpu kernels - to determine GPU block and grid sizes the iteration
                               region (i.e. field shape) has to be known. This can normally be inferred by the kernel
                               parameters - however in special cases like boundary conditions a manual specification
                               may be necessary.
    """
    assert isinstance(ghost_layers_to_include, str) or ghost_layers_to_include >= 0
    ast_params = kernel.parameters
    vec_info = ctx.get('cpu_vectorize_info', None)
    instruction_set = kernel.get_ast_attr('instruction_set')
    if vec_info:
        assume_inner_stride_one = vec_info['assume_inner_stride_one']
        assume_aligned = vec_info['assume_aligned']
        nontemporal = vec_info['nontemporal']
    else:
        assume_inner_stride_one = nontemporal = False
        assume_aligned = False

    cpu_openmp = ctx.get('cpu_openmp', False)
    kernel_ghost_layers = kernel.get_ast_attr('ghost_layers')

    ghost_layers_to_include = sp.sympify(ghost_layers_to_include)
    if kernel_ghost_layers is None:
        required_ghost_layers = 0
    else:
        # ghost layer info is ((x_gl_front, x_gl_end), (y_gl_front, y_gl_end).. )
        if isinstance(kernel_ghost_layers, int):
            required_ghost_layers = kernel_ghost_layers
        else:
            required_ghost_layers = max(max(kernel_ghost_layers))

    kernel_call_lines = []

    def get_cell_interval(field_object):
        if isinstance(cell_interval, str):
            return cell_interval
        elif isinstance(cell_interval, dict):
            return cell_interval[field_object]
        else:
            return None

    def get_start_coordinates(field_object):
        ci = get_cell_interval(field_object)
        if ci is None:
            return [-ghost_layers_to_include - required_ghost_layers] * field_object.spatial_dimensions
        else:
            assert ghost_layers_to_include == 0
            return [sp.Symbol(f"{ci}.{coord_name}Min()") - required_ghost_layers for coord_name in ('x', 'y', 'z')]

    def get_end_coordinates(field_object):
        ci = get_cell_interval(field_object)
        if ci is None:
            shape_names = ['xSize()', 'ySize()', 'zSize()'][:field_object.spatial_dimensions]
            offset = 2 * ghost_layers_to_include + 2 * required_ghost_layers
            return [f"int64_c({field_object.name}->{e}) + {offset}" for e in shape_names]
        else:
            assert ghost_layers_to_include == 0
            return [f"int64_c({ci}.{coord_name}Size()) + {2 * required_ghost_layers}"
                    for coord_name in ('x', 'y', 'z')]

    for param in ast_params:
        if param.is_field_parameter and FieldType.is_indexed(param.fields[0]):
            continue

        if param.is_field_pointer:
            field = param.fields[0]
            if field.field_type == FieldType.BUFFER:
                kernel_call_lines.append(f"{param.symbol.dtype} {param.symbol.name} = {param.field_name};")
            else:
                coordinates = get_start_coordinates(field)
                actual_gls = f"int_c({param.field_name}->nrOfGhostLayers())"
                coord_set = set(coordinates)
                coord_set = sorted(coord_set, key=lambda e: str(e))
                for c in coord_set:
                    kernel_call_lines.append(f"WALBERLA_ASSERT_GREATER_EQUAL({c}, -{actual_gls})")
                while len(coordinates) < 4:
                    coordinates.append(0)
                coordinates = tuple(coordinates)
                kernel_call_lines.append(f"{param.symbol.dtype} {param.symbol.name} = {param.field_name}->dataAt"
                                         f"({coordinates[0]}, {coordinates[1]}, {coordinates[2]}, {coordinates[3]});")
                if assume_inner_stride_one and field.index_dimensions > 0:
                    kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL({param.field_name}->layout(), field::fzyx)")
                if instruction_set and assume_aligned:
                    if nontemporal and cpu_openmp and 'cachelineZero' in instruction_set:
                        kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL((uintptr_t) {field.name}->dataAt(0, 0, 0, 0) %"
                                                 f"{instruction_set['cachelineSize']}, 0)")
                    else:
                        kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL((uintptr_t) {field.name}->dataAt(0, 0, 0, 0) %"
                                                 f"{instruction_set['bytes']}, 0)")
        elif param.is_field_stride:
            casted_stride = get_field_stride(param)
            type_str = param.symbol.dtype.c_name
            kernel_call_lines.append(f"const {type_str} {param.symbol.name} = {casted_stride};")
        elif param.is_field_shape:
            coord = param.symbol.coordinate
            field = param.fields[0]
            type_str = param.symbol.dtype.c_name
            shape = f"{type_str}({get_end_coordinates(field)[coord]})"
            assert coord < 3
            max_value = f"{field.name}->{('x', 'y', 'z')[coord]}SizeWithGhostLayer()"
            kernel_call_lines.append(f"WALBERLA_ASSERT_GREATER_EQUAL({max_value}, {shape})")
            kernel_call_lines.append(f"const {type_str} {param.symbol.name} = {shape};")
            if assume_inner_stride_one and field.index_dimensions > 0:
                kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL({field.name}->layout(), field::fzyx)")
            if instruction_set and assume_aligned:
                if nontemporal and cpu_openmp and 'cachelineZero' in instruction_set:
                    kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL((uintptr_t) {field.name}->dataAt(0, 0, 0, 0) %"
                                             f"{instruction_set['cachelineSize']}, 0)")
                else:
                    kernel_call_lines.append(f"WALBERLA_ASSERT_EQUAL((uintptr_t) {field.name}->dataAt(0, 0, 0, 0) %"
                                             f"{instruction_set['bytes']}, 0)")

    kernel_call_lines.append(kernel.generate_kernel_invocation_code(stream=stream,
                                                                    spatial_shape_symbols=spatial_shape_symbols))

    return "\n".join(kernel_call_lines)


@jinja2_context_decorator
def generate_function_collection_call(ctx, kernel, parameters_to_ignore=(),
                                      cell_interval=None, ghost_layers=None, use_field_ids=False):

    """Generates the function call to a pystencils kernel. It can be understood as a lightweight version of
       `generate_call`. Thus, it will only generate the parameters needed to call the kernel as a list of strings.

    Args:
        ctx: code generation context
        kernel: pystencils kernel
        parameters_to_ignore: In some cases not all parameters need to be printed. This is especially the case when
                              fixed parameters exist that are hardcoded in the jinja template.
        cell_interval: Defines the name (string) of a walberla CellInterval object in scope.
        ghost_layers: Defines the name (string) of a variable to define the number of used ghost_layers.
        use_field_ids: If set to true field names will be printed with the suffix `ID_`, to indicated that
                       a BlockDataID is passed.
    """

    target = translate_target(ctx['target'])
    is_gpu = target == Target.GPU

    parameters = []
    for param in kernel.parameters:
        if param.is_field_pointer and param.field_name not in parameters_to_ignore:
            if use_field_ids:
                parameters.append(f"{param.field_name}ID_")
            else:
                parameters.append(param.field_name)

    for param in kernel.parameters:
        if not param.is_field_parameter and param.symbol.name not in parameters_to_ignore:
            parameters.append(param.symbol.name)

    # TODO due to backward compatibility with high level interface spec
    for parameter in kernel.kernel_selection_tree.get_selection_parameter_list():
        if parameter.name not in parameters_to_ignore:
            parameters.append(parameter.name)

    if cell_interval:
        assert ghost_layers is None, "If a cell interval is specified ghost layers can not be specified"
        parameters.append(cell_interval)

    if ghost_layers:
        parameters.append(ghost_layers)

    if is_gpu and "gpuStream" not in parameters_to_ignore:
        parameters.append(f"gpuStream")

    return ", ".join(parameters)


def generate_swaps(kernel_info):
    """Generates code to swap main fields with temporary fields"""
    swaps = ""
    for src, dst in kernel_info.field_swaps:
        swaps += f"{src}->swapDataPointers({dst});\n"
    return swaps


def generate_timestep_advancements(kernel_info, advance=True):
    """Generates code to detect even or odd timestep"""
    if kernel_info.field_timestep:
        field_name = kernel_info.field_timestep["field_name"]
        advancement_function = kernel_info.field_timestep["function"]
        if advancement_function == "advanceTimestep" and advance is False:
            advancement_function = "getTimestepPlusOne"
        return f"uint8_t timestep = {field_name}->{advancement_function}();"
    return ""


def generate_constructor_initializer_list(kernel_infos, parameters_to_ignore=None, parameter_registration=None):
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]

    parameters_to_skip = []
    if parameters_to_ignore is not None:
        parameters_to_skip = [p for p in parameters_to_ignore]

    for kernel_info in kernel_infos:
        parameters_to_skip += kernel_info.temporary_fields

    parameter_initializer_list = []
    # First field pointer
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if param.is_field_pointer and param.field_name not in parameters_to_skip:
                parameter_initializer_list.append(f"{param.field_name}ID({param.field_name}ID_)")
                parameters_to_skip.append(param.field_name)

    # Then free parameters
    if parameter_registration is not None:
        parameters_to_skip.extend(parameter_registration.scaling_info)

    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if not param.is_field_parameter and param.symbol.name not in parameters_to_skip:
                parameter_initializer_list.append(f"{param.symbol.name}_({param.symbol.name})")
                parameters_to_skip.append(param.symbol.name)

    return ", ".join(parameter_initializer_list)


# TODO check varying_parameters
def generate_constructor_parameters(kernel_infos, parameters_to_ignore=None):
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]

    parameters_to_skip = []
    if parameters_to_ignore is not None:
        parameters_to_skip = [p for p in parameters_to_ignore]

    varying_parameters = []
    for kernel_info in kernel_infos:
        if hasattr(kernel_info, 'varying_parameters'):
            varying_parameters = kernel_info.varying_parameters
        varying_parameter_names = tuple(e[1] for e in varying_parameters)
        parameters_to_skip += kernel_info.temporary_fields + varying_parameter_names

    parameter_list = []
    # First field pointer
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if param.is_field_pointer and param.field_name not in parameters_to_skip:
                parameter_list.append(f"BlockDataID {param.field_name}ID_")
                parameters_to_skip.append(param.field_name)

    # Then free parameters
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if not param.is_field_parameter and param.symbol.name not in parameters_to_skip:
                parameter_list.append(f"{param.symbol.dtype} {param.symbol.name}")
                parameters_to_skip.append(param.symbol.name)

    varying_parameters = ["%s %s" % e for e in varying_parameters]
    return ", ".join(parameter_list + varying_parameters)


def generate_constructor_call_arguments(kernel_info, parameters_to_ignore=None):
    parameters_to_skip = []
    if parameters_to_ignore is not None:
        parameters_to_skip = [p for p in parameters_to_ignore]

    varying_parameters = []
    if hasattr(kernel_info, 'varying_parameters'):
        varying_parameters = kernel_info.varying_parameters
    varying_parameter_names = tuple(e[1] for e in varying_parameters)
    parameters_to_skip += kernel_info.temporary_fields + varying_parameter_names

    parameter_list = []
    for param in kernel_info.parameters:
        if param.is_field_pointer and param.field_name not in parameters_to_skip:
            parameter_list.append(f"{param.field_name}ID")
        elif not param.is_field_parameter and param.symbol.name not in parameters_to_skip:
            parameter_list.append(f'{param.symbol.name}_')
    varying_parameters = [f"{e}_" for e in varying_parameter_names]
    return ", ".join(parameter_list + varying_parameters)


@jinja2_context_decorator
def generate_members(ctx, kernel_infos, parameters_to_ignore=None, only_fields=False, parameter_registration=None):
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]

    if parameters_to_ignore is None:
        parameters_to_ignore = []

    params_to_skip = [p for p in parameters_to_ignore]

    fields = dict()
    for kernel_info in kernel_infos:
        for field in kernel_info.fields_accessed:
            fields[field.name] = field

    varying_parameters = []
    for kernel_info in kernel_infos:
        if hasattr(kernel_info, 'varying_parameters'):
            varying_parameters = kernel_info.varying_parameters
        varying_parameter_names = tuple(e[1] for e in varying_parameters)
        params_to_skip += kernel_info.temporary_fields
        params_to_skip += varying_parameter_names

    target = translate_target(ctx['target'])
    is_gpu = target == Target.GPU

    result = []
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if only_fields and not param.is_field_parameter:
                continue
            if param.is_field_pointer and param.field_name not in params_to_skip:
                result.append(f"BlockDataID {param.field_name}ID;")
                params_to_skip.append(param.field_name)

    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if only_fields and not param.is_field_parameter:
                continue
            if not param.is_field_parameter and param.symbol.name not in params_to_skip:
                if parameter_registration and param.symbol.name in parameter_registration.scaling_info:
                    result.append(f"std::vector<{param.symbol.dtype}> {param.symbol.name}Vector;")
                else:
                    result.append(f"{param.symbol.dtype} {param.symbol.name}_;")
                params_to_skip.append(param.symbol.name)

    for kernel_info in kernel_infos:
        for field_name in kernel_info.temporary_fields:
            f = fields[field_name]
            if field_name in parameters_to_ignore:
                continue
            parameters_to_ignore.append(field_name)
            assert field_name.endswith('_tmp')
            original_field_name = field_name[:-len('_tmp')]
            f_size = get_field_fsize(f)
            field_type = make_field_type(get_base_type(f.dtype), f_size, is_gpu)
            result.append(temporary_fieldMemberTemplate.format(type=field_type, original_field_name=original_field_name))

    for kernel_info in kernel_infos:
        if hasattr(kernel_info, 'varying_parameters'):
            result.extend(["%s %s_;" % e for e in kernel_info.varying_parameters])

    return "\n".join(result)


@jinja2_context_decorator
def generate_plain_parameter_list(ctx, kernel_info, cell_interval=None, ghost_layers=None, stream=None):
    fields = {f.name: f for f in kernel_info.fields_accessed}
    target = translate_target(ctx['target'])
    is_gpu = target == Target.GPU

    result = []
    for param in kernel_info.parameters:
        if not param.is_field_parameter:
            continue
        if param.is_field_pointer and param.field_name:
            f = fields[param.field_name]
            f_size = get_field_fsize(f)
            field_type = make_field_type(get_base_type(f.dtype), f_size, is_gpu)
            result.append(f"{field_type} * {param.field_name}")

    for param in kernel_info.parameters:
        if not param.is_field_parameter and param.symbol.name:
            result.append(f"{param.symbol.dtype} {param.symbol.name}")

    if hasattr(kernel_info, 'varying_parameters'):
        result.extend(["%s %s_;" % e for e in kernel_info.varying_parameters])

    # TODO due to backward compatibility with high level interface spec
    for parameter in kernel_info.kernel_selection_tree.get_selection_parameter_list():
        result.append(f"{parameter.dtype} {parameter.name}")

    if cell_interval:
        result.append(f"const CellInterval & {cell_interval}")

    if ghost_layers is not None:
        if type(ghost_layers) in (int, ):
            result.append(f"const cell_idx_t ghost_layers = {ghost_layers}")
        else:
            result.append(f"const cell_idx_t ghost_layers")

    if is_gpu:
        if stream is not None:
            result.append(f"gpuStream_t stream = {stream}")
        else:
            result.append(f"gpuStream_t stream")

    return ", ".join(result)


def generate_destructor(kernel_infos, class_name):
    temporary_fields = []
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]
    for kernel_info in kernel_infos:
        for tmp_field in kernel_info.temporary_fields:
            if tmp_field not in temporary_fields:
                temporary_fields.append(tmp_field)

    if not temporary_fields:
        return ""
    else:
        contents = ""
        for field_name in temporary_fields:
            contents += delete_loop.format(original_field_name=field_name[:-len('_tmp')])
        return temporary_constructor.format(contents=contents, class_name=class_name)


# IMPORTANT REMARK:
# This is specifically implemented for using generated kernels in the waLBerla's free surface LBM and is
# implemented in rather unflexible fashion. Therefore, it should not be extended and in the long-term, the free
# surface implementation should be refactored such that the general generated stream() is applicable.
@jinja2_context_decorator
def generate_field_type(ctx, kernel_info):
    fields = {f.name: f for f in kernel_info.fields_accessed}
    target = translate_target(ctx['target'])
    is_gpu = target == Target.GPU

    result = []
    for field_name in kernel_info.temporary_fields:
        f = fields[field_name]
        assert field_name.endswith('_tmp')
        original_field_name = field_name[:-len('_tmp')]
        f_size = get_field_fsize(f)
        field_type = make_field_type(get_base_type(f.dtype), f_size, is_gpu)
        result.append(temporary_fieldPointerTemplate.format(type=field_type, original_field_name=original_field_name))
    return "\n".join(result)


@jinja2_context_decorator
def nested_class_method_definition_prefix(ctx, nested_class_name):
    outer_class = ctx['class_name']
    if len(nested_class_name) == 0:
        return outer_class
    else:
        return f"{outer_class}::{nested_class_name}"


@jinja2_context_decorator
def generate_parameter_registration(ctx, kernel_infos, parameter_registration):
    if parameter_registration is None:
        return ""
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]

    params_to_skip = []
    result = []
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if not param.is_field_parameter and param.symbol.name not in params_to_skip:
                if param.symbol.name in parameter_registration.scaling_info:
                    result.append(standard_parameter_registration.format(dtype=param.symbol.dtype,
                                                                         name=param.symbol.name))
                    params_to_skip.append(param.symbol.name)

    return "\n".join(result)


@jinja2_context_decorator
def generate_constructor(ctx, kernel_infos, parameter_registration):
    if parameter_registration is None:
        return ""
    if not isinstance(kernel_infos, Iterable):
        kernel_infos = [kernel_infos]

    params_to_skip = []
    result = []
    for kernel_info in kernel_infos:
        for param in kernel_info.parameters:
            if not param.is_field_parameter and param.symbol.name not in params_to_skip:
                if param.symbol.name in parameter_registration.scaling_info:
                    name = param.symbol.name
                    dtype = param.symbol.dtype
                    result.append(standard_parameter_registration.format(dtype=dtype, name=name))
                    params_to_skip.append(name)

    return "\n".join(result)


def generate_list_of_expressions(expressions, prepend=''):
    if len(expressions) == 0:
        return ''
    return prepend + ", ".join(expressions)


def type_identifier_list(nested_arg_list):
    """
    Filters a nested list of strings and TypedSymbols and returns a comma-separated string.
    Strings are passed through as they are, but TypedSymbols are formatted as C-style
    'type identifier' strings, e.g. 'uint32_t ghost_layers'.
    """
    result = []

    def recursive_flatten(arg_list):
        for s in arg_list:
            if isinstance(s, str) and len(s) > 0:
                result.append(s)
            elif isinstance(s, TypedSymbol):
                result.append(f"{s.dtype} {s.name}")
            else:
                recursive_flatten(s)

    recursive_flatten(nested_arg_list)
    return ", ".join(result)


def identifier_list(nested_arg_list):
    """
    Filters a nested list of strings and TypedSymbols and returns a comma-separated string.
    Strings are passed through as they are, but TypedSymbols are replaced by their name.
    """
    result = []

    def recursive_flatten(arg_list):
        for s in arg_list:
            if isinstance(s, str):
                result.append(s)
            elif isinstance(s, TypedSymbol):
                result.append(s.name)
            else:
                recursive_flatten(s)

    recursive_flatten(nested_arg_list)
    return ", ".join(result)


def add_pystencils_filters_to_jinja_env(jinja_env):
    jinja_env.filters['generate_definition'] = generate_definition
    jinja_env.filters['generate_declaration'] = generate_declaration
    jinja_env.filters['generate_definitions'] = generate_definitions
    jinja_env.filters['generate_declarations'] = generate_declarations
    jinja_env.filters['generate_members'] = generate_members
    jinja_env.filters['generate_plain_parameter_list'] = generate_plain_parameter_list
    jinja_env.filters['generate_constructor_parameters'] = generate_constructor_parameters
    jinja_env.filters['generate_constructor_initializer_list'] = generate_constructor_initializer_list
    jinja_env.filters['generate_constructor_call_arguments'] = generate_constructor_call_arguments
    jinja_env.filters['generate_call'] = generate_call
    jinja_env.filters['generate_function_collection_call'] = generate_function_collection_call
    jinja_env.filters['generate_block_data_to_field_extraction'] = generate_block_data_to_field_extraction
    jinja_env.filters['generate_timestep_advancements'] = generate_timestep_advancements
    jinja_env.filters['generate_swaps'] = generate_swaps
    jinja_env.filters['generate_refs_for_kernel_parameters'] = generate_refs_for_kernel_parameters
    jinja_env.filters['generate_destructor'] = generate_destructor
    jinja_env.filters['generate_field_type'] = generate_field_type
    jinja_env.filters['nested_class_method_definition_prefix'] = nested_class_method_definition_prefix
    jinja_env.filters['generate_parameter_registration'] = generate_parameter_registration
    jinja_env.filters['generate_constructor'] = generate_constructor
    jinja_env.filters['type_identifier_list'] = type_identifier_list
    jinja_env.filters['identifier_list'] = identifier_list
    jinja_env.filters['list_of_expressions'] = generate_list_of_expressions
    jinja_env.filters['field_type'] = field_type
