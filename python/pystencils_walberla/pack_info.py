from collections import OrderedDict, defaultdict
from dataclasses import replace
from itertools import product
from typing import Dict, Optional, Sequence, Tuple

from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils import Assignment, AssignmentCollection, Field, FieldType, Target, create_kernel
from pystencils.backends.cbackend import get_headers
from pystencils.stencil import inverse_direction, offset_to_direction_string

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.kernel_info import KernelInfo
from pystencils_walberla.utility import config_from_context


def generate_pack_info_for_field(generation_context: CodeGenerationContext, class_name: str, field: Field,
                                 direction_subset: Optional[Tuple[Tuple[int, int, int]]] = None,
                                 operator=None, gl_to_inner=False,
                                 target=Target.CPU, data_type=None, cpu_openmp=False,
                                 **create_kernel_params):
    """Creates a pack info for a pystencils field assuming a pull-type stencil, packing all cell elements.

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        field: pystencils field for which to generate pack info
        direction_subset: optional sequence of directions for which values should be packed
                          otherwise a D3Q27 stencil is assumed
        operator: optional operator for, e.g., reduction pack infos
        gl_to_inner: communicates values from ghost layers of sender to interior of receiver
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """

    if not direction_subset:
        direction_subset = tuple((i, j, k) for i, j, k in product(*[(-1, 0, 1)] * 3))

    all_index_accesses = [field(*ind) for ind in product(*[range(s) for s in field.index_shape])]
    return generate_pack_info(generation_context, class_name, {direction_subset: all_index_accesses}, operator=operator,
                              gl_to_inner=gl_to_inner, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                              **create_kernel_params)


def generate_pack_info_from_kernel(generation_context: CodeGenerationContext, class_name: str,
                                   assignments: Sequence[Assignment], kind='pull', operator=None,
                                   target=Target.CPU, data_type=None, cpu_openmp=False,
                                   **create_kernel_params):
    """Generates a waLBerla GPU PackInfo from a (pull) kernel.

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        assignments: list of assignments from the compute kernel - generates PackInfo for "pull" part only
                     i.e. the kernel is expected to only write to the center
        kind: can either be pull or push
        operator: optional operator for, e.g., reduction pack infos
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
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
            for comm_dir in _comm_directions(comm_direction):
                spec[(comm_dir,)].add(fa.field.center(*fa.index))
    elif kind == 'push':
        for fa in writes:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            for comm_dir in _comm_directions(fa.offsets):
                spec[(comm_dir,)].add(fa)
    else:
        raise ValueError("Invalid 'kind' parameter")
    return generate_pack_info(generation_context, class_name, spec, operator=operator,
                              target=target, data_type=data_type, cpu_openmp=cpu_openmp, **create_kernel_params)


def generate_pack_info(generation_context: CodeGenerationContext, class_name: str,
                       directions_to_pack_terms: Dict[Tuple[Tuple], Sequence[Field.Access]],
                       namespace='pystencils', operator=None, gl_to_inner=False,
                       target=Target.CPU, data_type=None, cpu_openmp=False,
                       **create_kernel_params):
    """Generates a waLBerla GPU PackInfo

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name: name of the generated class
        directions_to_pack_terms: maps tuples of directions to read field accesses, specifying which values have to be
                                  packed for which direction
        namespace: inner namespace of the generated class
        operator: optional operator for, e.g., reduction pack infos
        gl_to_inner: communicates values from ghost layers of sender to interior of receiver
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    if cpu_openmp:
        raise ValueError("The packing kernels are already called inside an OpenMP parallel region. Thus "
                         "additionally parallelising each kernel is not supported.")
    items = [(e[0], sorted(e[1], key=lambda x: str(x))) for e in directions_to_pack_terms.items()]
    items = sorted(items, key=lambda e: e[0])
    directions_to_pack_terms = OrderedDict(items)

    config = config_from_context(generation_context, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 **create_kernel_params)

    config_zero_gl = config_from_context(generation_context, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                         ghost_layers=0, **create_kernel_params)

    # Vectorisation of the pack info is not implemented.
    config = replace(config, cpu_vectorize_info=None)
    config_zero_gl = replace(config_zero_gl, cpu_vectorize_info=None)

    config = replace(config, allow_double_writes=True)
    config_zero_gl = replace(config_zero_gl, allow_double_writes=True)

    template_name = "CpuPackInfo.tmpl" if config.target == Target.CPU else 'GpuPackInfo.tmpl'

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
        err_detail = "\n".join(f" - {f.name} [{f.dtype}]" for f in fields_accessed)
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
        pack_ast = create_kernel(pack_assignments, config=config_zero_gl)
        pack_ast.function_name = 'pack_{}'.format("_".join(direction_strings))
        if operator is None:
            unpack_assignments = [Assignment(term, buffer(i)) for i, term in enumerate(terms)]
        else:
            unpack_assignments = [Assignment(term, operator(term, buffer(i))) for i, term in enumerate(terms)]
        unpack_ast = create_kernel(unpack_assignments, config=config_zero_gl)
        unpack_ast.function_name = 'unpack_{}'.format("_".join(direction_strings))

        pack_kernels[direction_strings] = KernelInfo(pack_ast)
        unpack_kernels[direction_strings] = KernelInfo(unpack_ast)
        elements_per_cell[direction_strings] = len(terms)
    fused_kernel = create_kernel([Assignment(buffer.center, t) for t in all_accesses], config=config)

    jinja_context = {
        'class_name': class_name,
        'pack_kernels': pack_kernels,
        'unpack_kernels': unpack_kernels,
        'fused_kernel': KernelInfo(fused_kernel),
        'elements_per_cell': elements_per_cell,
        'headers': get_headers(fused_kernel),
        'target': config.target.name.lower(),
        'dtype': dtype,
        'field_name': field_names.pop(),
        'namespace': namespace,
        'gl_to_inner': gl_to_inner,
    }
    env = Environment(loader=PackageLoader('pystencils_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    header = env.get_template(template_name + ".h").render(**jinja_context)
    source = env.get_template(template_name + ".cpp").render(**jinja_context)

    source_extension = "cu" if target == Target.GPU and generation_context.cuda else "cpp"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)


def generate_mpidtype_info_from_kernel(ctx: CodeGenerationContext, class_name: str,
                                       assignments: Sequence[Assignment], kind='pull', namespace='pystencils'):
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
            for comm_dir in _comm_directions(comm_direction):
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
            for comm_dir in _comm_directions(fa.offsets):
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
    ctx.write_file(f"{class_name}.h", header)


# ---------------------------------- Internal --------------------------------------------------------------------------

def _comm_directions(direction):
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
