from os import path
from functools import reduce
from typing import Union, Dict, DefaultDict
import warnings

from pystencils import CreateKernelConfig, Target
from pystencils.backends.simd_instruction_sets import get_supported_instruction_sets
from pystencils.boundaries.createindexlist import boundary_index_array_coordinate_names, direction_member_name
from pystencils.typing import BasicType, create_type, get_base_type

from lbmpy import LBStencil

from pystencils_walberla.cmake_integration import CodeGenerationContext

HEADER_EXTENSIONS = {'.h', '.hpp'}


def generate_info_header(ctx: CodeGenerationContext,
                         filename: str,
                         stencil_typedefs: dict = None,
                         field_typedefs: dict = None,
                         additional_headers: set = None,
                         headers_to_ignore: set = None,
                         additional_typedefs: dict = None,
                         additional_code: str = ""):
    """Generates an info header, consolidating required information about the generated code.
    The info header #includes all generated header files, and is thus the only header the
    application needs to #include. It can also contain aliases for waLBerla stencil types and
    instances of the GhostLayerField template.

    Args:
        ctx: Code Generation Context
        filename: Name of the generated info header file
        stencil_typedefs: dict mapping type names to stencil names or tuples
        field_typedefs: dict mapping type names to pystencils `Field` instances
        additional_headers: additional header files to be included
        headers_to_ignore: headers which should not be included
        additional_typedefs: dict mapping aliases to types
        additional_code: additional code which gets appended on the file
    """
    stencil_typedefs = stencil_typedefs if stencil_typedefs is not None else dict()
    field_typedefs = field_typedefs if field_typedefs is not None else dict()
    additional_typedefs = additional_typedefs if additional_typedefs is not None else dict()

    additional_headers = additional_headers if additional_headers is not None else set()
    headers_to_ignore = headers_to_ignore if headers_to_ignore is not None else set()

    headers_in_ctx = set(_filter_headers(ctx.files_written))

    stencil_headers, stencil_typedefs = _stencil_inclusion_code(stencil_typedefs)
    field_headers, field_typedefs = _field_inclusion_code(field_typedefs)

    headers_to_include = stencil_headers | field_headers | headers_in_ctx | additional_headers
    headers_to_include = sorted(headers_to_include - headers_to_ignore)
    headers_to_include = [f'#include "{header}"' for header in headers_to_include]

    typedefs = {**stencil_typedefs, **field_typedefs, **additional_typedefs}
    typedefs = [f"using {alias} = {typename};" for alias, typename in typedefs.items()]

    lines = "#pragma once\n"

    lines += '\n'.join(headers_to_include + [''] + typedefs) + '\n'

    if path.splitext(filename)[1] not in HEADER_EXTENSIONS:
        filename += '.h'

    ctx.write_file(filename, lines + additional_code)


def get_vectorize_instruction_set(ctx: CodeGenerationContext):
    """returns a list of supported vector instruction sets. If waLBerla is not build with
       `WALBERLA_OPTIMIZE_FOR_LOCALHOST` `None` is returned.

    Args:
        ctx: Code Generation Context
    """

    if ctx.optimize_for_localhost:
        supported_instruction_sets = get_supported_instruction_sets()
        if supported_instruction_sets:
            return supported_instruction_sets[-1]
        else:  # if cpuinfo package is not installed
            warnings.warn("Could not obtain supported vectorization instruction sets - defaulting to sse. "
                          "This problem can probably be fixed by installing py-cpuinfo. This package can "
                          "gather the needed hardware information.")
            return 'sse'
    else:
        return None


def config_from_context(ctx: CodeGenerationContext, target: Target = Target.CPU,
                        data_type: Union[type, str, DefaultDict[str, BasicType], Dict[str, BasicType]] = None,
                        cpu_openmp: Union[bool, int] = None, cpu_vectorize_info: Dict = None,
                        **kwargs) -> CreateKernelConfig:
    """Creates a :class: `pystencils.config.CreateKernelConfig` from the code generation context. By default,
       all arguments are determined by the generation context. This means for example if `DWALBERLA_BUILD_WITH_GPU_SUPPORT` is
       `True` the kernel will be generated for GPU using either CUDA or HIP.

    Args:
        ctx: Code Generation Context
        target: All targets are defined in :class:`pystencils.enums.Target`
        data_type: Data type used for all untyped symbols (i.e. non-fields), can also be a dict from symbol name to
                   type. If specified as a dict ideally a defaultdict is used to define a default value for symbols
                   not listed in the dict. If a plain dict is provided it will be transformed into a defaultdict
                   internally. The default value will then be specified via type collation then.
        cpu_openmp: `True` or number of threads for OpenMP parallelization, `False` for no OpenMP.
                     If set to `True`, the maximum number of available threads will be chosen.
        cpu_vectorize_info: A dictionary with keys, 'vector_instruction_set', 'assume_aligned' and 'nontemporal'
                            for documentation of these parameters see vectorize function. Example:
                            '{'instruction_set': 'avx512', 'assume_aligned': True, 'nontemporal':True}'
        kwargs: keyword arguments that can be taken by :class: `pystencils.config.CreateKernelConfig`
    """

    if target == Target.GPU and not ctx.gpu:
        raise ValueError("can not generate gpu code if waLBerla is not build with GPU support. Please use "
                         "-DWALBERLA_BUILD_WITH_CUDA=1 or -DWALBERLA_BUILD_WITH_HIP=1 for configuring cmake")

    default_dtype = "float64" if ctx.double_accuracy else "float32"
    if data_type is None:
        data_type = default_dtype

    if cpu_openmp and not ctx.openmp:
        warnings.warn("Code is generated with OpenMP pragmas but waLBerla is not build with OpenMP. "
                      "The compilation might not work due to wrong compiler flags. "
                      "Please use -DWALBERLA_BUILD_WITH_OPENMP=1 for configuring cmake")

    if cpu_openmp is None:
        cpu_openmp = ctx.openmp

    if cpu_vectorize_info is None:
        cpu_vectorize_info = {}

    default_vec_is = get_vectorize_instruction_set(ctx)

    cpu_vectorize_info['instruction_set'] = cpu_vectorize_info.get('instruction_set', default_vec_is)
    cpu_vectorize_info['assume_inner_stride_one'] = cpu_vectorize_info.get('assume_inner_stride_one', True)
    cpu_vectorize_info['assume_aligned'] = cpu_vectorize_info.get('assume_aligned', False)
    cpu_vectorize_info['nontemporal'] = cpu_vectorize_info.get('nontemporal', False)
    cpu_vectorize_info['assume_sufficient_line_padding'] = cpu_vectorize_info.get('assume_sufficient_line_padding',
                                                                                  False)

    config = CreateKernelConfig(target=target, data_type=data_type, default_number_float=data_type,
                                cpu_openmp=cpu_openmp, cpu_vectorize_info=cpu_vectorize_info,
                                **kwargs)

    return config


def merge_sorted_lists(lx, ly, sort_key=lambda x: x, identity_check_key=None):
    if identity_check_key is None:
        identity_check_key = sort_key
    nx = len(lx)
    ny = len(ly)

    def recursive_merge(lx_intern, ly_intern, ix_intern, iy_intern):
        if ix_intern == nx:
            return ly_intern[iy_intern:]
        if iy_intern == ny:
            return lx_intern[ix_intern:]
        x = lx_intern[ix_intern]
        y = ly_intern[iy_intern]
        skx = sort_key(x)
        sky = sort_key(y)
        if skx == sky:
            if identity_check_key(x) == identity_check_key(y):
                return [x] + recursive_merge(lx_intern, ly_intern, ix_intern + 1, iy_intern + 1)
            else:
                raise ValueError(f'Elements <{x}> and <{y}> with equal sort key where not identical!')
        elif skx < sky:
            return [x] + recursive_merge(lx_intern, ly_intern, ix_intern + 1, iy_intern)
        else:
            return [y] + recursive_merge(lx_intern, ly_intern, ix_intern, iy_intern + 1)
    return recursive_merge(lx, ly, 0, 0)


def merge_lists_of_symbols(lists):
    def merger(lx, ly):
        return merge_sorted_lists(lx, ly, sort_key=lambda x: x.symbol.name, identity_check_key=lambda x: x.symbol)
    return reduce(merger, lists)


def struct_from_numpy_dtype(struct_name, numpy_dtype):
    result = f"struct {struct_name} {{ \n"

    equality_compare = []
    constructor_params = []
    constructor_initializer_list = []
    for name, (sub_type, offset) in numpy_dtype.fields.items():
        pystencils_type = create_type(sub_type)
        result += f"    {pystencils_type} {name};\n"
        if name in boundary_index_array_coordinate_names or name == direction_member_name:
            constructor_params.append(f"{pystencils_type} {name}_")
            constructor_initializer_list.append(f"{name}({name}_)")
        else:
            constructor_initializer_list.append(f"{name}()")
        if pystencils_type.is_float():
            equality_compare.append(f"floatIsEqual({name}, o.{name})")
        else:
            equality_compare.append(f"{name} == o.{name}")

    result += "    %s(%s) : %s {}\n" % \
              (struct_name, ", ".join(constructor_params), ", ".join(constructor_initializer_list))
    result += "    bool operator==(const %s & o) const {\n        return %s;\n    }\n" % \
              (struct_name, " && ".join(equality_compare))
    result += "};\n"
    return result


#   ------------------------------------- INTERNAL -------------------------------------------------------------


def _filter_headers(filepaths):
    for p in filepaths:
        if path.splitext(p)[1] in HEADER_EXTENSIONS:
            yield path.split(p)[1]


def _stencil_inclusion_code(stencil_typedefs):
    headers = set()
    typedefs = dict()
    for typename, stencil in stencil_typedefs.items():
        if isinstance(stencil, tuple):
            dim = len(stencil[0])
            q = len(stencil)
            stencil = f"D{dim}Q{q}"
        elif isinstance(stencil, LBStencil):
            stencil = stencil.name
        elif not isinstance(stencil, str):
            raise ValueError(f'Invalid stencil: Do not know what to make of {stencil}')

        headers.add(f"stencil/{stencil}.h")
        typedefs[typename] = f"walberla::stencil::{stencil}"

    return headers, typedefs


def _field_inclusion_code(field_typedefs):
    headers = set()
    typedefs = dict()
    for typename, field in field_typedefs.items():
        f_size = field.values_per_cell()
        dtype = get_base_type(field.dtype)
        headers.add("field/GhostLayerField.h")
        typedefs[typename] = f"walberla::field::GhostLayerField<{dtype}, {f_size}>"

    return headers, typedefs
