import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model

with CodeGeneration() as ctx:
    omega_shear, omega_bulk = sp.symbols("omega_shear, omega_bulk")
    temperature = sp.symbols("temperature")
    force_field, vel_field = ps.fields("force(3), velocity(3): [3D]", layout='fzyx')

    collision_rule = create_lb_collision_rule(
        stencil='D3Q19', compressible=True, fluctuating={
            'temperature' : temperature,
            'block_offsets' : 'walberla',
        },
        method='mrt3', relaxation_rates=[omega_shear, omega_bulk],
        force_model='guo', force=force_field.center_vector,
        optimization={'cse_global': True}
    )

    generate_lattice_model(ctx, 'FluctuatingMRT_LatticeModel', collision_rule)
