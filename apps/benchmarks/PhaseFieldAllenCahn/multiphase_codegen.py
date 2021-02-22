from pystencils import fields, TypedSymbol
from pystencils.simp import sympy_cse
from pystencils import AssignmentCollection

from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.stencils import get_stencil

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel

from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, \
    hydrodynamic_force, get_collision_assignments_hydro

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


force_h = [f / 3 for f in interface_tracking_force(C, stencil_phase, W)]
force_model_h = MultiphaseForceModel(force=force_h)

force_g = hydrodynamic_force(g, C, method_hydro, relaxation_time, density_liquid, density_gas, kappa, beta, body_force)

h_tmp_symbol_list = [h_tmp.center(i) for i, _ in enumerate(stencil_phase)]
sum_h = np.sum(h_tmp_symbol_list[:])

####################
# LBM UPDATE RULES #
####################

method_phase.set_force_model(force_model_h)

phase_field_LB_step = create_lb_update_rule(lb_method=method_phase,
                                            velocity_input=u,
                                            compressible=True,
                                            optimization={"symbolic_field": h,
                                                          "symbolic_temporary_field": h_tmp},
                                            kernel_type='stream_pull_collide')


phase_field_LB_step.set_main_assignments_from_dict({**phase_field_LB_step.main_assignments_dict, **{C.center: sum_h}})

phase_field_LB_step = AssignmentCollection(main_assignments=phase_field_LB_step.main_assignments,
                                           subexpressions=phase_field_LB_step.subexpressions)
phase_field_LB_step = sympy_cse(phase_field_LB_step)

# ---------------------------------------------------------------------------------------------------------

hydro_LB_step = get_collision_assignments_hydro(lb_method=method_hydro,
                                                density=density,
                                                velocity_input=u,
                                                force=force_g,
                                                sub_iterations=1,
                                                symbolic_fields={"symbolic_field": g,
                                                                 "symbolic_temporary_field": g_tmp},
                                                kernel_type='collide_stream_push')

# streaming of the hydrodynamic distribution
stream_hydro = create_lb_update_rule(stencil=stencil_hydro,
                                     optimization={"symbolic_field": g,
                                                   "symbolic_temporary_field": g_tmp},
                                     kernel_type='stream_pull_only')

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
                       field_swaps=[(h, h_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec)

        generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                       field_swaps=[(g, g_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec)

        # communication
        generate_pack_info_from_kernel(ctx, 'PackInfo_phase_field_distributions',
                                       phase_field_LB_step.main_assignments, target='cpu')
        generate_pack_info_from_kernel(ctx, 'PackInfo_phase_field',
                                       hydro_LB_step.all_assignments, target='cpu', kind='pull')
        generate_pack_info_from_kernel(ctx, 'PackInfo_velocity_based_distributions',
                                       hydro_LB_step.all_assignments, target='cpu', kind='push')

        ctx.write_file("GenDefines.h", info_header)

    if ctx.cuda:
        generate_sweep(ctx, 'initialize_phase_field_distributions',
                       h_updates, target='gpu')
        generate_sweep(ctx, 'initialize_velocity_based_distributions',
                       g_updates, target='gpu')

        generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                       field_swaps=[(h, h_tmp)],
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
        generate_pack_info_from_kernel(ctx, 'PackInfo_phase_field_distributions',
                                       phase_field_LB_step.main_assignments, target='gpu')
        generate_pack_info_from_kernel(ctx, 'PackInfo_phase_field',
                                       hydro_LB_step.all_assignments, target='gpu', kind='pull')
        generate_pack_info_from_kernel(ctx, 'PackInfo_velocity_based_distributions',
                                       hydro_LB_step.all_assignments, target='gpu', kind='push')

        ctx.write_file("GenDefines.h", info_header)

print("finished code generation successfully")
