import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule, create_mrt_orthogonal, force_model_from_string
from lbmpy.moments import is_bulk_moment, is_shear_moment, get_order
from lbmpy.stencils import get_stencil
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model

with CodeGeneration() as ctx:
    omega_shear = sp.symbols("omega_shear")
    temperature = sp.symbols("temperature")
    force_field, vel_field = ps.fields("force(3), velocity(3): [3D]", layout='fzyx')

    def rr_getter(moment_group):
        is_shear = [is_shear_moment(m, 3) for m in moment_group]
        is_bulk = [is_bulk_moment(m, 3) for m in moment_group]
        order = [get_order(m) for m in moment_group]
        assert min(order) == max(order)
        order = order[0]

        if order < 2:
            return 0
        elif any(is_bulk):
            assert all(is_bulk)
            return sp.Symbol("omega_bulk")
        elif any(is_shear):
            assert all(is_shear)
            return sp.Symbol("omega_shear")
        elif order % 2 == 0:
            assert order > 2
            return sp.Symbol("omega_even")
        else:
            return sp.Symbol("omega_odd")

    method = create_mrt_orthogonal(
        stencil=get_stencil('D3Q19'),
        compressible=True,
        weighted=True,
        relaxation_rate_getter=rr_getter,
        force_model=force_model_from_string('schiller', force_field.center_vector)
    )
    collision_rule = create_lb_collision_rule(
        method,
        fluctuating={
            'temperature': temperature,
            'block_offsets': 'walberla',
            'rng_node': ps.rng.PhiloxTwoDoubles if ctx.double_accuracy else ps.rng.PhiloxFourFloats,
        },
        optimization={'cse_global': True}
    )

    params = {}
    if ctx.optimize_for_localhost:
        params['cpu_vectorize_info'] = {'assume_inner_stride_one': True, 'assume_aligned': True}
    generate_lattice_model(ctx, 'FluctuatingMRT_LatticeModel', collision_rule, field_layout='fzyx', **params)
