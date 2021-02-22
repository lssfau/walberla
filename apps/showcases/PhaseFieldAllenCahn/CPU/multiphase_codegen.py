from pystencils import fields
from pystencils.simp import sympy_cse
from pystencils import AssignmentCollection

from lbmpy.boundaries import NoSlip
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.stencils import get_stencil

from pystencils_walberla import CodeGeneration, generate_sweep
from lbmpy_walberla import generate_boundary

from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, \
    hydrodynamic_force, get_collision_assignments_hydro

from lbmpy.phasefield_allen_cahn.force_model import MultiphaseForceModel

import numpy as np
import sympy as sp
import pystencils as ps

stencil_phase = get_stencil("D3Q19")
stencil_hydro = get_stencil("D3Q27")
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
W = 5
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

method_phase = create_lb_method(stencil=stencil_phase, method='srt',
                                relaxation_rate=relaxation_rate_allen_cahn, compressible=True)

method_hydro = create_lb_method(stencil=stencil_hydro, method="mrt", weighted=True,
                                relaxation_rates=[relaxation_rate, 1, 1, 1, 1, 1])

# create the kernels for the initialization of the g and h field
h_updates = initializer_kernel_phase_field_lb(h, C, u, method_phase, W, fd_stencil=get_stencil("D3Q27"))
g_updates = initializer_kernel_hydro_lb(g, u, method_hydro)

force_h = [f / 3 for f in interface_tracking_force(C, stencil_phase, W, fd_stencil=get_stencil("D3Q27"))]
force_model_h = MultiphaseForceModel(force=force_h)

force_g = hydrodynamic_force(g, C, method_hydro, relaxation_time, density_liquid, density_gas, kappa, beta, body_force,
                             fd_stencil=get_stencil("D3Q27"))

####################
# LBM UPDATE RULES #
####################

h_tmp_symbol_list = [h_tmp.center(i) for i, _ in enumerate(stencil_phase)]
sum_h = np.sum(h_tmp_symbol_list[:])

method_phase.set_force_model(force_model_h)

phase_field_LB_step = create_lb_update_rule(lb_method=method_phase,
                                            velocity_input=u,
                                            compressible=True,
                                            optimization={"symbolic_field": h,
                                                          "symbolic_temporary_field": h_tmp},
                                            kernel_type='stream_pull_collide')

phase_field_LB_step.set_main_assignments_from_dict({**phase_field_LB_step.main_assignments_dict, **{C.center: sum_h}})
subexp = [ps.Assignment(a.lhs, float(a.rhs)) if a.rhs == 0 else a for a in phase_field_LB_step.subexpressions]
phase_field_LB_step = AssignmentCollection(main_assignments=phase_field_LB_step.main_assignments,
                                           subexpressions=subexp)
phase_field_LB_step = sympy_cse(phase_field_LB_step)

# ---------------------------------------------------------------------------------------------------------
hydro_LB_step = get_collision_assignments_hydro(lb_method=method_hydro,
                                                density=density,
                                                velocity_input=u,
                                                force=force_g,
                                                sub_iterations=2,
                                                symbolic_fields={"symbolic_field": g,
                                                                 "symbolic_temporary_field": g_tmp},
                                                kernel_type='collide_only')

hydro_LB_step.set_sub_expressions_from_dict({**{relaxation_rate: relaxation_rate_cutoff},
                                             **hydro_LB_step.subexpressions_dict})

stream_hydro = create_lb_update_rule(stencil=stencil_hydro,
                                     optimization={"symbolic_field": g,
                                                   "symbolic_temporary_field": g_tmp},
                                     kernel_type='stream_pull_only')

###################
# GENERATE SWEEPS #
###################
cpu_vec = {'instruction_set': 'sse', 'assume_inner_stride_one': True, 'nontemporal': True}

vp = [('int32_t', 'cudaBlockSize0'),
      ('int32_t', 'cudaBlockSize1')]

info_header = f"""
#include "stencil/D3Q{q_phase}.h"\nusing Stencil_phase_T = walberla::stencil::D3Q{q_phase};
#include "stencil/D3Q{q_hydro}.h"\nusing Stencil_hydro_T = walberla::stencil::D3Q{q_hydro};
"""

with CodeGeneration() as ctx:
    if not ctx.optimize_for_localhost:
        cpu_vec = {'instruction_set': None}

    generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates, target='cpu')
    generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates, target='cpu')

    generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                   field_swaps=[(h, h_tmp)],
                   inner_outer_split=True,
                   cpu_vectorize_info=cpu_vec, target='cpu')
    generate_boundary(ctx, 'phase_field_LB_NoSlip', NoSlip(), method_phase, target='cpu')

    generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                   inner_outer_split=True,
                   cpu_vectorize_info=cpu_vec, target='cpu')
    generate_boundary(ctx, 'hydro_LB_NoSlip', NoSlip(), method_hydro, target='cpu')

    generate_sweep(ctx, 'stream_hydro', stream_hydro,
                   field_swaps=[(g, g_tmp)],
                   inner_outer_split=True,
                   cpu_vectorize_info=cpu_vec, target='cpu')

    ctx.write_file("GenDefines.h", info_header)

    # TODO: generate communication (PackInfo)

print("finished code generation successfully")
