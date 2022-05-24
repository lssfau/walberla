import sympy as sp

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule

from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model

#   ========================
#      General Parameters
#   ========================

stencil = LBStencil(Stencil.D2Q9)
omega = sp.Symbol('omega')
layout = 'fzyx'

#   Optimizations for the LBM Method
lbm_opt = LBMOptimisation(cse_global=True, field_layout=layout)

#   ===========================
#      SRT Method Definition
#   ===========================

lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=omega,
                       zero_centered=False, delta_equilibrium=False)

srt_collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
srt_update_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:
    # generation of the lattice model ...
    generate_lattice_model(ctx, "SRTLatticeModel", srt_collision_rule, field_layout=layout)
    # ... and generation of the pack information to be used for the MPI communication
    generate_pack_info_from_kernel(ctx, "SRTPackInfo", srt_update_rule)
