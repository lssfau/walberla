import sympy as sp

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy import LBMConfig, LBMOptimisation, Method, Stencil
from lbmpy.stencils import LBStencil

from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model

# general parameters
stencil = LBStencil(Stencil.D3Q19)
omega = sp.Symbol('omega')
layout = 'fzyx'

# optimizations to be used by the code generator
optimizations = {'cse_global': True, 'field_layout': layout}

# method definition
method_params = {'stencil': stencil,
                 'method': 'srt',
                 'relaxation_rate': omega,
                 'compressible': True}

lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=omega, compressible=True,
                       zero_centered=False, delta_equilibrium=False,
                       streaming_pattern='pull')
lbm_opt = LBMOptimisation(cse_global=True, field_layout=layout)

collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

with CodeGeneration() as ctx:
    generate_lattice_model(ctx, "GeneratedLatticeModel", collision_rule, field_layout=layout)
