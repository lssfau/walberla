from collections import defaultdict
from lbmpy import LBStencil, Stencil
from lbmpy.advanced_streaming.utility import Timestep, get_accessor, get_timesteps
from lbmpy.advanced_streaming.communication import _extend_dir

from pystencils import Assignment, Field, Target
from pystencils.stencil import inverse_direction

from pystencils_walberla.pack_info import _comm_directions, generate_pack_info


def generate_lb_pack_info(generation_context,
                          class_name_prefix: str,
                          stencil,
                          pdf_field,
                          streaming_pattern='pull',
                          lb_collision_rule=None,
                          always_generate_separate_classes=False,
                          namespace='lbm',
                          target=Target.CPU,
                          data_type=None,
                          cpu_openmp=False,
                          **create_kernel_params):
    """Generates waLBerla MPI PackInfos for an LBM kernel, based on a given method
    and streaming pattern. For in-place streaming patterns, two PackInfos are generated;
    one for the even and another for the odd time steps.

    Args:
        generation_context: see documentation of `generate_sweep`
        class_name_prefix: Prefix of the desired class name which will be extended with
                           'Even' or 'Odd' for in-place kernels
        stencil: Instance of class LBStencil
        pdf_field: pdf field for which the pack info is created
        streaming_pattern: The employed streaming pattern.
        lb_collision_rule: Optional. The collision rule defining the lattice boltzmann kernel, as returned
                           by `create_lb_collision_rule`. If specified, it will be scanned for non-local
                           accesses to other fields other than the PDF fields (as might be required for
                           computing gradients in coupled simulations), whose communication will then
                           be included in the PackInfo.
        always_generate_separate_classes: If True, generate a pair of Even/Odd PackInfos even for a two-
                                          fields kernel (i.e. the pull/push patterns). Otherwise, for two-fields
                                          kernels, only one PackInfo class will be generated without a
                                          suffix to its name.
        namespace: inner namespace of the generated class
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
        **create_kernel_params: remaining keyword arguments are passed to `pystencils.create_kernel`
    """
    timesteps = [Timestep.EVEN, Timestep.ODD] \
        if always_generate_separate_classes \
        else get_timesteps(streaming_pattern)

    common_spec = defaultdict(set)

    if lb_collision_rule is not None:
        assignments = lb_collision_rule.all_assignments
        reads = set()
        for a in assignments:
            if not isinstance(a, Assignment):
                continue
            reads.update(a.rhs.atoms(Field.Access))
        for fa in reads:
            assert all(abs(e) <= 1 for e in fa.offsets)
            if all(offset == 0 for offset in fa.offsets):
                continue
            comm_direction = inverse_direction(fa.offsets)
            for comm_dir in _comm_directions(comm_direction):
                common_spec[(comm_dir,)].add(fa.field.center(*fa.index))

    full_stencil = LBStencil(Stencil.D3Q27) if stencil.D == 3 else LBStencil(Stencil.D2Q9)

    for t in timesteps:
        spec = common_spec.copy()
        write_accesses = get_accessor(streaming_pattern, t).write(pdf_field, stencil)
        for comm_dir in full_stencil:
            if all(d == 0 for d in comm_dir):
                continue

            for streaming_dir in set(_extend_dir(comm_dir)) & set(stencil):
                d = stencil.index(streaming_dir)
                fa = write_accesses[d]
                spec[(comm_dir,)].add(fa)

        if t == Timestep.EVEN:
            class_name_suffix = 'Even'
        elif t == Timestep.ODD:
            class_name_suffix = 'Odd'
        else:
            class_name_suffix = ''

        class_name = class_name_prefix + class_name_suffix
        generate_pack_info(generation_context, class_name, spec, namespace=namespace,
                           target=target, data_type=data_type, cpu_openmp=cpu_openmp, **create_kernel_params)
