import sympy as sp
import pystencils as ps

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule
from lbmpy.boundaries import NoSlip

from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model, generate_boundary

#   =====================
#   Common Parameters
#   =====================

STENCIL = 'D2Q9'
OMEGA = sp.Symbol('omega')
LAYOUT = 'zyxf'

#   Optimization
OPT = {'target': 'cpu', 'cse_global': True, 'field_layout': LAYOUT}

#   Velocity Output Field
vel_field = ps.fields("velocity(2): [2D]", layout=LAYOUT)
OUTPUT = {'velocity': vel_field}

#   ==================================================================================================
#   Method Definitions, including collision rules, boundary handling and communication pack infos.
#   ==================================================================================================


#   SRT Method
def build_srt_method(ctx):
    srt_params = {'stencil': STENCIL,
                  'method': 'srt',
                  'relaxation_rate': OMEGA}

    srt_collision_rule = create_lb_collision_rule(optimization=OPT, output=OUTPUT, **srt_params)
    generate_lattice_model(ctx, "SRTLatticeModel", srt_collision_rule, field_layout=LAYOUT)

    srt_update_rule = create_lb_update_rule(collision_rule=srt_collision_rule, optimization=OPT)
    generate_pack_info_from_kernel(ctx, "SRTPackInfo", srt_update_rule)

    generate_boundary(ctx, "SRTNoSlip", NoSlip(), srt_collision_rule.method)


#   Cumulant MRT Method
def build_cumulant_method(ctx):
    mrt_cumulant_params = {'stencil': STENCIL,
                           'method': 'mrt_raw',
                           'relaxation_rates': [0, 0, 0, OMEGA, OMEGA, OMEGA, 1, 1, 1],
                           'cumulant': True,
                           'compressible': True}   # Incompressible cumulants not yet supported!

    mrt_cumulant_collision_rule = create_lb_collision_rule(optimization=OPT, output=OUTPUT, **mrt_cumulant_params)
    generate_lattice_model(ctx, "CumulantMRTLatticeModel", mrt_cumulant_collision_rule, field_layout=LAYOUT)

    mrt_cumulant_update_rule = create_lb_update_rule(collision_rule=mrt_cumulant_collision_rule, optimization=OPT)
    generate_pack_info_from_kernel(ctx, "CumulantMRTPackInfo", mrt_cumulant_update_rule)

    generate_boundary(ctx, "CumulantMRTNoSlip", NoSlip(), mrt_cumulant_collision_rule.method)


#   Orthogonal MRT Method with entropy constraint
def build_entropic_method(ctx):
    omega_f = sp.Symbol('omega_f')

    mrt_entropic_params = {'stencil': STENCIL,
                           'method': 'mrt',
                           'relaxation_rates': [OMEGA, OMEGA, omega_f, omega_f],
                           'entropic': True,
                           'compressible': True}    # Entropic models only implemented with pdfs centered around 1

    mrt_entropic_collision_rule = create_lb_collision_rule(optimization=OPT, output=OUTPUT, **mrt_entropic_params)
    generate_lattice_model(ctx, "EntropicMRTLatticeModel", mrt_entropic_collision_rule, field_layout=LAYOUT)

    mrt_entropic_update_rule = create_lb_update_rule(collision_rule=mrt_entropic_collision_rule, optimization=OPT)
    generate_pack_info_from_kernel(ctx, "EntropicMRTPackInfo", mrt_entropic_update_rule)

    generate_boundary(ctx, "EntropicMRTNoSlip", NoSlip(), mrt_entropic_collision_rule.method)

#   ================================
#   Main Function
#   ================================


if __name__ == "__main__":
    #   All code generation must happen within one context, and all files specified in
    #   CMakeLists.txt must be generated therein.
    with CodeGeneration() as ctx:
        build_srt_method(ctx)
        build_cumulant_method(ctx)
        build_entropic_method(ctx)
