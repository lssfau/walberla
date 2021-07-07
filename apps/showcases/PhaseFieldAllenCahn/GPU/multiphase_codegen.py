from lbmpy.phasefield_allen_cahn.contact_angle import ContactAngle
from pystencils import fields, TypedSymbol
from pystencils.simp import sympy_cse
from pystencils import Assignment
from pystencils.astnodes import Block, Conditional

from lbmpy.boundaries import NoSlip
from lbmpy.creationfunctions import create_lb_method
from lbmpy.stencils import get_stencil

import pystencils_walberla
from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_for_field
from lbmpy_walberla import generate_boundary, generate_lb_pack_info

from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb,\
    initializer_kernel_hydro_lb, interface_tracking_force, hydrodynamic_force, get_collision_assignments_hydro,\
    get_collision_assignments_phase

from lbmpy.phasefield_allen_cahn.force_model import MultiphaseForceModel

import numpy as np
import sympy as sp

stencil_phase_name = "D3Q27"
stencil_hydro_name = "D3Q27"

contact_angle_in_degrees = 22

stencil_phase = get_stencil(stencil_phase_name)
stencil_hydro = get_stencil(stencil_hydro_name)
q_phase = len(stencil_phase)
q_hydro = len(stencil_hydro)

assert (len(stencil_phase[0]) == len(stencil_hydro[0]))
dimensions = len(stencil_hydro[0])

########################
# PARAMETER DEFINITION #
########################

density_liquid = sp.Symbol("rho_H")
density_gas = sp.Symbol("rho_L")

surface_tension = sp.Symbol("sigma")
mobility = sp.Symbol("mobility")

gravitational_acceleration = sp.Symbol("gravity")

relaxation_time_liquid = sp.Symbol("tau_H")
relaxation_time_gas = sp.Symbol("tau_L")

# phase-field parameter
drho3 = (density_liquid - density_gas) / 3
# interface thickness
W = sp.Symbol("interface_thickness")
# coefficients related to surface tension
beta = 12.0 * (surface_tension / W)
kappa = 1.5 * surface_tension * W

########################
# FIELDS #
########################

# velocity field
u = fields(f"vel_field({dimensions}): [{dimensions}D]", layout='fzyx')
# phase-field
C = fields(f"phase_field: [{dimensions}D]", layout='fzyx')
C_tmp = fields(f"phase_field_tmp: [{dimensions}D]", layout='fzyx')

flag = fields(f"flag_field: uint8[{dimensions}D]", layout='fzyx')
# phase-field distribution functions
h = fields(f"lb_phase_field({q_phase}): [{dimensions}D]", layout='fzyx')
h_tmp = fields(f"lb_phase_field_tmp({q_phase}): [{dimensions}D]", layout='fzyx')
# hydrodynamic distribution functions
g = fields(f"lb_velocity_field({q_hydro}): [{dimensions}D]", layout='fzyx')
g_tmp = fields(f"lb_velocity_field_tmp({q_hydro}): [{dimensions}D]", layout='fzyx')

########################################
# RELAXATION RATES AND EXTERNAL FORCES #
########################################

# relaxation rate for interface tracking LB step
relaxation_rate_allen_cahn = 1.0 / (0.5 + (3.0 * mobility))
# calculate the relaxation rate for hydrodynamic LB step
density = density_gas + C.center * (density_liquid - density_gas)
# force acting on all phases of the model
body_force = np.array([0, gravitational_acceleration * density, 0])
# calculation of the relaxation time via viscosities
# viscosity = viscosity_gas * viscosity_liquid + C.center\
#             * (density_liquid*viscosity_liquid - viscosity_liquid*viscosity_gas)
# relaxation_time = 3 * viscosity / density + 0.5

relaxation_time = 0.5 + relaxation_time_gas + C.center * (relaxation_time_liquid - relaxation_time_gas)
# calculate the relaxation time if the phase-field has values over one
# relaxation_rate = 1.0 / relaxation_time
relaxation_rate = sp.Symbol("s8")
relaxation_rate_cutoff = sp.Piecewise((1 / (0.5 + relaxation_time_liquid), C.center > 0.999),   # True value
                                      (1 / relaxation_time, True))                              # Else value

###############
# LBM METHODS #
###############

# method_phase = create_lb_method(stencil=stencil_phase, method="mrt", compressible=True, weighted=True,
#                                 relaxation_rates=[1, 1.5, 1, 1.5, 1, 1.5])
method_phase = create_lb_method(stencil=stencil_phase, method="mrt", compressible=True, weighted=True,
                                relaxation_rates=[1, 1, 1, 1, 1, 1])

method_phase.set_conserved_moments_relaxation_rate(relaxation_rate_allen_cahn)

method_hydro = create_lb_method(stencil=stencil_hydro, method="mrt", weighted=True,
                                relaxation_rates=[relaxation_rate, 1, 1, 1, 1, 1])


# create the kernels for the initialization of the g and h field
h_updates = initializer_kernel_phase_field_lb(h, C, u, method_phase, W, fd_stencil=get_stencil("D3Q27"))
g_updates = initializer_kernel_hydro_lb(g, u, method_hydro)

force_h = [f / 3 for f in interface_tracking_force(C, stencil_phase, W, fd_stencil=get_stencil("D3Q27"))]
force_model_h = MultiphaseForceModel(force=force_h)

force_g = hydrodynamic_force(g, C, method_hydro, relaxation_time, density_liquid, density_gas, kappa, beta, body_force,
                             fd_stencil=get_stencil("D3Q27"))

force_model_g = MultiphaseForceModel(force=force_g, rho=density)

####################
# LBM UPDATE RULES #
####################

phase_field_LB_step = get_collision_assignments_phase(lb_method=method_phase,
                                                      velocity_input=u,
                                                      output={'density': C_tmp},
                                                      force_model=force_model_h,
                                                      symbolic_fields={"symbolic_field": h,
                                                                       "symbolic_temporary_field": h_tmp},
                                                      kernel_type='stream_pull_collide')

phase_field_LB_step = sympy_cse(phase_field_LB_step)

phase_field_LB_step = [Conditional(sp.Eq(flag.center(), 2),
                                   Block(phase_field_LB_step),
                                   Block([Assignment(C_tmp.center, C.center)]))]
# ---------------------------------------------------------------------------------------------------------
hydro_LB_step = get_collision_assignments_hydro(lb_method=method_hydro,
                                                density=density,
                                                velocity_input=u,
                                                force_model=force_model_g,
                                                sub_iterations=2,
                                                symbolic_fields={"symbolic_field": g,
                                                                 "symbolic_temporary_field": g_tmp},
                                                kernel_type='collide_stream_push')

hydro_LB_step.set_sub_expressions_from_dict({**{relaxation_rate: relaxation_rate_cutoff},
                                             **hydro_LB_step.subexpressions_dict})

hydro_LB_step = [Conditional(sp.Eq(flag.center(), 2), Block(hydro_LB_step))]

contact_angle = ContactAngle(contact_angle_in_degrees, W)


###################
# GENERATE SWEEPS #
###################

vp = [('int32_t', 'cudaBlockSize0'),
      ('int32_t', 'cudaBlockSize1'),
      ('int32_t', 'cudaBlockSize2')]

sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                    TypedSymbol("cudaBlockSize1", np.int32),
                    TypedSymbol("cudaBlockSize2", np.int32))

sweep_params = {'block_size': sweep_block_size}

info_header = f"""
using namespace walberla;
#include "stencil/D3Q{q_phase}.h"\nusing Stencil_phase_T = walberla::stencil::D3Q{q_phase};
#include "stencil/D3Q{q_hydro}.h"\nusing Stencil_hydro_T = walberla::stencil::D3Q{q_hydro};
using PdfField_phase_T = GhostLayerField<real_t, {q_phase}>;
using PdfField_hydro_T = GhostLayerField<real_t, {q_hydro}>;
using VelocityField_T = GhostLayerField<real_t, {dimensions}>;
using PhaseField_T = GhostLayerField<real_t, 1>;
#ifndef UTIL_H
#define UTIL_H
const char * stencil_phase_name = "{stencil_phase_name}";
const char * stencil_hydro_name = "{stencil_hydro_name}";
#endif
"""

with CodeGeneration() as ctx:
    generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates, target='gpu')
    generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates, target='gpu')

    generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                   field_swaps=[(h, h_tmp), (C, C_tmp)],
                   inner_outer_split=True,
                   target='gpu',
                   gpu_indexing_params=sweep_params,
                   varying_parameters=vp)
    generate_boundary(ctx, 'phase_field_LB_NoSlip', NoSlip(), method_phase, target='gpu', streaming_pattern='pull')

    generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                   field_swaps=[(g, g_tmp)],
                   inner_outer_split=True,
                   target='gpu',
                   gpu_indexing_params=sweep_params,
                   varying_parameters=vp)
    generate_boundary(ctx, 'hydro_LB_NoSlip', NoSlip(), method_hydro, target='gpu', streaming_pattern='push')

    # communication

    generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                          streaming_pattern='pull', target='gpu')

    generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                          streaming_pattern='push', target='gpu')

    generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target='gpu')

    pystencils_walberla.boundary.generate_boundary(ctx, 'ContactAngle', contact_angle,
                                                   C.name, stencil_hydro, index_shape=[], target='gpu')

    ctx.write_file("GenDefines.h", info_header)

print("finished code generation successfully")
