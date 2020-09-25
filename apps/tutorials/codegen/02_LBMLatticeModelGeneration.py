import sympy as sp

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule

from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model

#   ========================
#      General Parameters
#   ========================

stencil = 'D2Q9'
omega = sp.Symbol('omega')
layout = 'fzyx'

#   Optimizations to be used by the code generator
optimizations = {'cse_global': True, 'field_layout': layout}

#   ===========================
#      SRT Method Definition
#   ===========================

srt_params = {'stencil': stencil,
              'method': 'srt',
              'relaxation_rate': omega}

srt_collision_rule = create_lb_collision_rule(optimization=optimizations, **srt_params)
srt_update_rule = create_lb_update_rule(collision_rule=srt_collision_rule, optimization=optimizations)

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:
    # generation of the lattice model ...
    generate_lattice_model(ctx, "SRTLatticeModel", srt_collision_rule, field_layout=layout)
    # ... and generation of the pack information to be used for the MPI communication
    generate_pack_info_from_kernel(ctx, "SRTPackInfo", srt_update_rule)
