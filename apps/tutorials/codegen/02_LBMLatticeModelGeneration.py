import sympy as sp

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule

from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model

#   ========================
#      General Parameters
#   ========================

STENCIL = 'D2Q9'
OMEGA = sp.Symbol('omega')
LAYOUT = 'fzyx'

#   Optimization
OPT = {'target': 'cpu', 'cse_global': True, 'field_layout': LAYOUT}

#   ===========================
#      SRT Method Definition
#   ===========================

srt_params = {'stencil': STENCIL,
              'method': 'srt',
              'relaxation_rate': OMEGA}

srt_collision_rule = create_lb_collision_rule(optimization=OPT, **srt_params)
srt_update_rule = create_lb_update_rule(collision_rule=srt_collision_rule, optimization=OPT)

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:
    generate_lattice_model(ctx, "SRTLatticeModel", srt_collision_rule, field_layout=LAYOUT)
    generate_pack_info_from_kernel(ctx, "SRTPackInfo", srt_update_rule)
