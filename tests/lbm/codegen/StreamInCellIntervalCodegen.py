import sympy as sp

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule
from lbmpy import LBMConfig

from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model

# general parameters
stencil = 'D3Q19'
omega = sp.Symbol('omega')
layout = 'fzyx'

# optimizations to be used by the code generator
optimizations = {'cse_global': True, 'field_layout': layout}

# method definition
method_params = {'stencil': stencil,
                 'method': 'srt',
                 'relaxation_rate': omega,
                 'compressible': True}

lbm_config = LBMConfig(streaming_pattern='pull')

collision_rule = create_lb_collision_rule(optimization=optimizations, **method_params)
update_rule = create_lb_update_rule(collision_rule=collision_rule, lbm_config=lbm_config, optimization=optimizations)

with CodeGeneration() as ctx:
    generate_lattice_model(ctx, "GeneratedLatticeModel", collision_rule, field_layout=layout)
