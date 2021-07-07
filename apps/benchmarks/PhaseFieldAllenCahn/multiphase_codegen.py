from pystencils import fields, TypedSymbol
from pystencils.simp import sympy_cse
from pystencils import AssignmentCollection

from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.stencils import get_stencil

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_for_field
from lbmpy_walberla import generate_lb_pack_info

from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, \
    hydrodynamic_force, get_collision_assignments_hydro, get_collision_assignments_phase

from lbmpy.phasefield_allen_cahn.force_model import MultiphaseForceModel

import numpy as np

stencil_phase = get_stencil("D3Q15", "walberla")
stencil_hydro = get_stencil("D3Q27", "walberla")
q_phase = len(stencil_phase)
q_hydro = len(stencil_hydro)

maxwellian_moments = False
assert (len(stencil_phase[0]) == len(stencil_hydro[0]))
dimensions = len(stencil_hydro[0])

########################
# PARAMETER DEFINITION #
########################

density_liquid = 1.0
density_gas = 0.001
surface_tension = 0.0001
M = 0.02

# phase-field parameter
drho3 = (density_liquid - density_gas) / 3
# interface thickness
W = 5
# coefficient related to surface tension
beta = 12.0 * (surface_tension / W)
# coefficient related to surface tension
kappa = 1.5 * surface_tension * W
# relaxation rate allen cahn (h)
w_c = 1.0 / (0.5 + (3.0 * M))

########################
# FIELDS #
########################

# velocity field
u = fields(f"vel_field({dimensions}): [{dimensions}D]", layout='fzyx')
# phase-field
C = fields(f"phase_field: [{dimensions}D]", layout='fzyx')
C_tmp = fields(f"phase_field_tmp: [{dimensions}D]", layout='fzyx')

# phase-field distribution functions
h = fields(f"lb_phase_field({q_phase}): [{dimensions}D]", layout='fzyx')
h_tmp = fields(f"lb_phase_field_tmp({q_phase}): [{dimensions}D]", layout='fzyx')
# hydrodynamic distribution functions
g = fields(f"lb_velocity_field({q_hydro}): [{dimensions}D]", layout='fzyx')
g_tmp = fields(f"lb_velocity_field_tmp({q_hydro}): [{dimensions}D]", layout='fzyx')

########################################
# RELAXATION RATES AND EXTERNAL FORCES #
########################################

# calculate the relaxation rate for the hydro lb as well as the body force
density = density_gas + C.center * (density_liquid - density_gas)

# force acting on all phases of the model
body_force = np.array([0, 1e-6, 0])

relaxation_time = 0.03 + 0.5
relaxation_rate = 1.0 / relaxation_time

###############
# LBM METHODS #
###############

method_phase = create_lb_method(stencil=stencil_phase, method='srt', relaxation_rate=w_c, compressible=True)

method_hydro = create_lb_method(stencil=stencil_hydro, method="mrt", weighted=True,
                                relaxation_rates=[relaxation_rate, 1, 1, 1, 1, 1],
                                maxwellian_moments=maxwellian_moments)

# create the kernels for the initialization of the g and h field
h_updates = initializer_kernel_phase_field_lb(h, C, u, method_phase, W)
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

# ---------------------------------------------------------------------------------------------------------

hydro_LB_step = get_collision_assignments_hydro(lb_method=method_hydro,
                                                density=density,
                                                velocity_input=u,
                                                force_model=force_model_g,
                                                sub_iterations=2,
                                                symbolic_fields={"symbolic_field": g,
                                                                 "symbolic_temporary_field": g_tmp},
                                                kernel_type='collide_stream_push')

###################
# GENERATE SWEEPS #
###################

cpu_vec = {'assume_inner_stride_one': True, 'nontemporal': True}

vp = [('int32_t', 'cudaBlockSize0'),
      ('int32_t', 'cudaBlockSize1')]

sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                    TypedSymbol("cudaBlockSize1", np.int32),
                    1)
sweep_params = {'block_size': sweep_block_size}

info_header = f"""
#include "stencil/D3Q{q_phase}.h"\nusing Stencil_phase_T = walberla::stencil::D3Q{q_phase};
#include "stencil/D3Q{q_hydro}.h"\nusing Stencil_hydro_T = walberla::stencil::D3Q{q_hydro};
"""

with CodeGeneration() as ctx:
    if not ctx.cuda:
        if not ctx.optimize_for_localhost:
            cpu_vec = {'instruction_set': None}

        generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates)
        generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates)

        generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                       field_swaps=[(h, h_tmp), (C, C_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec)

        generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                       field_swaps=[(g, g_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec)

        # communication
        generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                              streaming_pattern='pull', target='cpu')

        generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                              streaming_pattern='push', target='cpu')

        generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target='cpu')

        ctx.write_file("GenDefines.h", info_header)

    if ctx.cuda:
        generate_sweep(ctx, 'initialize_phase_field_distributions',
                       h_updates, target='gpu')
        generate_sweep(ctx, 'initialize_velocity_based_distributions',
                       g_updates, target='gpu')

        generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                       field_swaps=[(h, h_tmp), (C, C_tmp)],
                       inner_outer_split=True,
                       target='gpu',
                       gpu_indexing_params=sweep_params,
                       varying_parameters=vp)

        generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                       field_swaps=[(g, g_tmp)],
                       inner_outer_split=True,
                       target='gpu',
                       gpu_indexing_params=sweep_params,
                       varying_parameters=vp)
        # communication
        generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                              streaming_pattern='pull', target='gpu')

        generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                              streaming_pattern='push', target='gpu')

        generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target='gpu')

        ctx.write_file("GenDefines.h", info_header)

print("finished code generation successfully")
