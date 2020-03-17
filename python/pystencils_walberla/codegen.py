from collections import OrderedDict, defaultdict
from itertools import product
from typing import Dict, Optional, Sequence, Tuple

from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils import (
    Assignment, AssignmentCollection, Field, FieldType, create_kernel, create_staggered_kernel)
from pystencils.astnodes import KernelFunction
from pystencils.backends.cbackend import get_headers
from pystencils.backends.simd_instruction_sets import get_supported_instruction_sets
from pystencils.stencil import inverse_direction, offset_to_direction_string
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env

__all__ = ['generate_sweep', 'generate_pack_info', 'generate_pack_info_for_field', 'generate_pack_info_from_kernel',
           'generate_mpidtype_info_from_kernel', 'default_create_kernel_parameters', 'KernelInfo']


def generate_sweep(generation_context, class_name, assignments,
                   namespace='pystencils', field_swaps=(), staggered=False, varying_parameters=(),
                   inner_outer_split=False,
                   **create_kernel_params):
    """Generates a waLBerla sweep from a pystencils representation.

    The constructor of the C++ sweep class expects all kernel parameters (fields and parameters) in alphabetical order.
    Fields have to passed using BlockDataID's pointing to walberla fields

    Args:
        generation_context: build system context filled with information from waLBerla's CMake. The context for example
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
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    create_kernel_params = default_create_kernel_parameters(generation_context, create_kernel_params)

    if not generation_context.cuda and create_kernel_params['target'] == 'gpu':
        return

    if isinstance(assignments, KernelFunction):
        ast = assignments
        create_kernel_params['target'] = ast.target
    elif not staggered:
        ast = create_kernel(assignments, **create_kernel_params)
    else:
        ast = create_staggered_kernel(assignments, **create_kernel_params)

    def to_name(f):
        return f.name if isinstance(f, Field) else f

    field_swaps = tuple((to_name(e[0]), to_name(e[1])) for e in field_swaps)
    temporary_fields = tuple(e[1] for e in field_swaps)

    ast.function_name = class_name.lower()

    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    if inner_outer_split is False:
        jinja_context = {
            'kernel': KernelInfo(ast, temporary_fields, field_swaps, varying_parameters),
            'namespace': namespace,
            'class_name': class_name,
            'target': create_kernel_params.get("target", "cpu"),
            'headers': get_headers(ast),
        }
        header = env.get_template("Sweep.tmpl.h").render(**jinja_context)
        source = env.get_template("Sweep.tmpl.cpp").render(**jinja_context)
    else:
        main_kernel_info = KernelInfo(ast, temporary_fields, field_swaps, varying_parameters)
        representative_field = {p.field_name for p in main_kernel_info.parameters if p.is_field_parameter}
        representative_field = sorted(representative_field)[0]

        jinja_context = {
            'kernel': main_kernel_info,
            'namespace': namespace,
            'class_name': class_name,
            'target': create_kernel_params.get("target", "cpu"),
            'field': representative_field,
            'headers': get_headers(ast),
        }
        header = env.get_template("SweepInnerOuter.tmpl.h").render(**jinja_context)
        source = env.get_template("SweepInnerOuter.tmpl.cpp").render(**jinja_context)

    source_extension = "cpp" if create_kernel_params.get("target", "cpu") == "cpu" else "cu"
    generation_context.write_file("{}.h".format(class_name), header)
    generation_context.write_file("{}.{}".format(class_name, source_extension), source)


def generate_pack_info_for_field(generation_context, class_name: str, field: Field,
                                 direction_subset: Optional[Tuple[Tuple[int, int, int]]] = None,
                                 **create_kernel_params):
    """Creates a pack info for a pystencils field assuming a pull-type stencil, packing all cell elements.

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        field: pystencils field for which to generate pack info
        direction_subset: optional sequence of directions for which values should be packed
                          otherwise a D3Q27 stencil is assumed
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    if not direction_subset:
        direction_subset = tuple((i, j, k) for i, j, k in product(*[(-1, 0, 1)] * 3))

    all_index_accesses = [field(*ind) for ind in product(*[range(s) for s in field.index_shape])]
    return generate_pack_info(generation_context, class_name, {direction_subset: all_index_accesses},
                              **create_kernel_params)


def generate_pack_info_from_kernel(generation_context, class_name: str, assignments: Sequence[Assignment],
                                   kind='pull', **create_kernel_params):
    """Generates a waLBerla GPU PackInfo from a (pull) kernel.

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        assignments: list of assignments from the compute kernel - generates PackInfo for "pull" part only
                     i.e. the kernel is expected to only write to the center
        kind:                      
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    assert kind in ('push', 'pull')
    reads = set()
    writes = set()

    if isinstance(assignments, AssignmentCollection):
        assignments = assignments.all_assignments

    for a in assignments:
        if not isinstance(a, Assignment):
            continue
        reads.update(a.rhs.atoms(Field.Access))
        writes.update(a.lhs.atoms(Field.Access))
    spec = defaultdict(set)
    if kind == 'pull':
        for fa in reads:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            comm_direction = inverse_direction(fa.offsets)
            for comm_dir in comm_directions(comm_direction):
                spec[(comm_dir,)].add(fa.field.center(*fa.index))
    elif kind == 'push':
        for fa in writes:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            for comm_dir in comm_directions(fa.offsets):
                spec[(comm_dir,)].add(fa)
    else:
        raise ValueError("Invalid 'kind' parameter")
    return generate_pack_info(generation_context, class_name, spec, **create_kernel_params)


def generate_pack_info(generation_context, class_name: str,
                       directions_to_pack_terms: Dict[Tuple[Tuple], Sequence[Field.Access]],
                       namespace='pystencils',
                       **create_kernel_params):
    """Generates a waLBerla GPU PackInfo

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        directions_to_pack_terms: maps tuples of directions to read field accesses, specifying which values have to be
                                  packed for which direction
        namespace: inner namespace of the generated class
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    items = [(e[0], sorted(e[1], key=lambda x: str(x))) for e in directions_to_pack_terms.items()]
    items = sorted(items, key=lambda e: e[0])
    directions_to_pack_terms = OrderedDict(items)

    create_kernel_params = default_create_kernel_parameters(generation_context, create_kernel_params)
    target = create_kernel_params.get('target', 'cpu')

    template_name = "CpuPackInfo.tmpl" if target == 'cpu' else 'GpuPackInfo.tmpl'

    fields_accessed = set()
    for terms in directions_to_pack_terms.values():
        for term in terms:
            assert isinstance(term, Field.Access)  # and all(e == 0 for e in term.offsets)
            fields_accessed.add(term)

    field_names = {fa.field.name for fa in fields_accessed}

    data_types = {fa.field.dtype for fa in fields_accessed}
    if len(data_types) == 0:
        raise ValueError("No fields to pack!")
    if len(data_types) != 1:
        err_detail = "\n".join(" - {} [{}]".format(f.name, f.dtype) for f in fields_accessed)
        raise NotImplementedError("Fields of different data types are used - this is not supported.\n" + err_detail)
    dtype = data_types.pop()

    pack_kernels = OrderedDict()
    unpack_kernels = OrderedDict()
    all_accesses = set()
    elements_per_cell = OrderedDict()
    for direction_set, terms in directions_to_pack_terms.items():
        for d in direction_set:
            if not all(abs(i) <= 1 for i in d):
                raise NotImplementedError("Only first neighborhood supported")

        buffer = Field.create_generic('buffer', spatial_dimensions=1, field_type=FieldType.BUFFER,
                                      dtype=dtype.numpy_dtype, index_shape=(len(terms),))

        direction_strings = tuple(offset_to_direction_string(d) for d in direction_set)
        all_accesses.update(terms)

        pack_assignments = [Assignment(buffer(i), term) for i, term in enumerate(terms)]
        pack_ast = create_kernel(pack_assignments, **create_kernel_params, ghost_layers=0)
        pack_ast.function_name = 'pack_{}'.format("_".join(direction_strings))
        unpack_assignments = [Assignment(term, buffer(i)) for i, term in enumerate(terms)]
        unpack_ast = create_kernel(unpack_assignments, **create_kernel_params, ghost_layers=0)
        unpack_ast.function_name = 'unpack_{}'.format("_".join(direction_strings))

        pack_kernels[direction_strings] = KernelInfo(pack_ast)
        unpack_kernels[direction_strings] = KernelInfo(unpack_ast)
        elements_per_cell[direction_strings] = len(terms)

    fused_kernel = create_kernel([Assignment(buffer.center, t) for t in all_accesses], **create_kernel_params)

    jinja_context = {
        'class_name': class_name,
        'pack_kernels': pack_kernels,
        'unpack_kernels': unpack_kernels,
        'fused_kernel': KernelInfo(fused_kernel),
        'elements_per_cell': elements_per_cell,
        'headers': get_headers(fused_kernel),
        'target': target,
        'dtype': dtype,
        'field_name': field_names.pop(),
        'namespace': namespace,
    }
    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    header = env.get_template(template_name + ".h").render(**jinja_context)
    source = env.get_template(template_name + ".cpp").render(**jinja_context)

    source_extension = "cpp" if target == "cpu" else "cu"
    generation_context.write_file("{}.h".format(class_name), header)
    generation_context.write_file("{}.{}".format(class_name, source_extension), source)


def generate_mpidtype_info_from_kernel(generation_context, class_name: str,
                                       assignments: Sequence[Assignment], kind='pull', namespace='pystencils', ):
    assert kind in ('push', 'pull')
    reads = set()
    writes = set()

    if isinstance(assignments, AssignmentCollection):
        assignments = assignments.all_assignments

    for a in assignments:
        if not isinstance(a, Assignment):
            continue
        reads.update(a.rhs.atoms(Field.Access))
        writes.update(a.lhs.atoms(Field.Access))

    spec = defaultdict(set)
    if kind == 'pull':
        read_fields = set(fa.field for fa in reads)
        assert len(read_fields) == 1, "Only scenarios where one fields neighbors are accessed"
        field = read_fields.pop()
        for fa in reads:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            comm_direction = inverse_direction(fa.offsets)
            for comm_dir in comm_directions(comm_direction):
                assert len(fa.index) == 1, "Supports only fields with a single index dimension"
                spec[(offset_to_direction_string(comm_dir),)].add(fa.index[0])
    elif kind == 'push':
        written_fields = set(fa.field for fa in writes)
        assert len(written_fields) == 1, "Only scenarios where one fields neighbors are accessed"
        field = written_fields.pop()

        for fa in writes:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            for comm_dir in comm_directions(fa.offsets):
                assert len(fa.index) == 1, "Supports only fields with a single index dimension"
                spec[(offset_to_direction_string(comm_dir),)].add(fa.index[0])
    else:
        raise ValueError("Invalid 'kind' parameter")

    jinja_context = {
        'class_name': class_name,
        'namespace': namespace,
        'kind': kind,
        'field_name': field.name,
        'f_size': field.index_shape[0],
        'spec': spec,
    }
    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    header = env.get_template("MpiDtypeInfo.tmpl.h").render(**jinja_context)
    generation_context.write_file("{}.h".format(class_name), header)


# ---------------------------------- Internal --------------------------------------------------------------------------


class KernelInfo:
    def __init__(self, ast, temporary_fields=(), field_swaps=(), varying_parameters=()):
        self.ast = ast
        self.temporary_fields = tuple(temporary_fields)
        self.field_swaps = tuple(field_swaps)
        self.varying_parameters = tuple(varying_parameters)
        self.parameters = ast.get_parameters()  # cache parameters here


def default_create_kernel_parameters(generation_context, params):
    default_dtype = "float64" if generation_context.double_accuracy else 'float32'

    if generation_context.optimize_for_localhost:
        supported_instruction_sets = get_supported_instruction_sets()
        if supported_instruction_sets:
            default_vec_is = get_supported_instruction_sets()[-1]
        else:  # if cpuinfo package is not installed
            default_vec_is = 'sse'
    else:
        default_vec_is = None

    params['target'] = params.get('target', 'cpu')
    params['data_type'] = params.get('data_type', default_dtype)
    params['cpu_openmp'] = params.get('cpu_openmp', generation_context.openmp)
    params['cpu_vectorize_info'] = params.get('cpu_vectorize_info', {})

    vec = params['cpu_vectorize_info']
    vec['instruction_set'] = vec.get('instruction_set', default_vec_is)
    vec['assume_inner_stride_one'] = True
    vec['assume_aligned'] = vec.get('assume_aligned', False)
    vec['nontemporal'] = vec.get('nontemporal', False)
    return params


def comm_directions(direction):
    if all(e == 0 for e in direction):
        yield direction
    binary_numbers_list = binary_numbers(len(direction))
    for comm_direction in binary_numbers_list:
        for i in range(len(direction)):
            if direction[i] == 0:
                comm_direction[i] = 0
            if direction[i] == -1 and comm_direction[i] == 1:
                comm_direction[i] = -1
        if not all(e == 0 for e in comm_direction):
            yield tuple(comm_direction)


def binary_numbers(n):
    result = list()
    for i in range(1 << n):
        binary_number = bin(i)[2:]
        binary_number = '0' * (n - len(binary_number)) + binary_number
        result.append((list(map(int, binary_number))))
    return result
