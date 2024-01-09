import numpy as np
import sympy as sp

from pystencils import fields, Target
from pystencils.typing import TypedSymbol
from pystencils.simp import sympy_cse

from lbmpy import LBMConfig, LBStencil, Method, Stencil
from lbmpy.creationfunctions import LBMOptimisation, create_lb_method, create_lb_update_rule
from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, hydrodynamic_force, add_interface_tracking_force, \
    add_hydrodynamic_force, hydrodynamic_force_assignments
from lbmpy.phasefield_allen_cahn.parameter_calculation import AllenCahnParameters

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_for_field, generate_info_header
from lbmpy_walberla import generate_lb_pack_info

with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"

    stencil_phase = LBStencil(Stencil.D3Q15)
    stencil_hydro = LBStencil(Stencil.D3Q27)
    assert (stencil_phase.D == stencil_hydro.D)

    ########################
    # PARAMETER DEFINITION #
    ########################
    # In the codegneration skript we only need the symbolic parameters
    parameters = AllenCahnParameters(density_heavy=1.0, density_light=0.001,
                                     dynamic_viscosity_heavy=0.03, dynamic_viscosity_light=0.03,
                                     surface_tension=0.0001, mobility=0.02, interface_thickness=5,
                                     gravitational_acceleration=1e-6)

    ########################
    # FIELDS #
    ########################

    # velocity field
    u = fields(f"vel_field({stencil_hydro.D}): {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    # phase-field
    C = fields(f"phase_field: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    C_tmp = fields(f"phase_field_tmp: {field_type}[{stencil_hydro.D}D]", layout='fzyx')

    # phase-field distribution functions
    h = fields(f"lb_phase_field({stencil_phase.Q}): {field_type}[{stencil_phase.D}D]", layout='fzyx')
    h_tmp = fields(f"lb_phase_field_tmp({stencil_phase.Q}): {field_type}[{stencil_phase.D}D]", layout='fzyx')
    # hydrodynamic distribution functions
    g = fields(f"lb_velocity_field({stencil_hydro.Q}): {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    g_tmp = fields(f"lb_velocity_field_tmp({stencil_hydro.Q}): {field_type}[{stencil_hydro.D}D]", layout='fzyx')

    ########################################
    # RELAXATION RATES AND EXTERNAL FORCES #
    ########################################

    # relaxation rate for interface tracking LB step
    relaxation_rate_allen_cahn = 1.0 / (0.5 + (3.0 * parameters.symbolic_mobility))
    # calculate the relaxation rate for hydrodynamic LB step
    rho_L = parameters.symbolic_density_light
    rho_H = parameters.symbolic_density_heavy
    density = rho_L + C.center * (rho_H - rho_L)
    # force acting on all phases of the model
    body_force = [0, 0, 0]
    body_force[1] = parameters.symbolic_gravitational_acceleration * density

    omega = parameters.omega(C)

    ###############
    # LBM METHODS #
    ###############

    rates = [0.0]
    rates += [relaxation_rate_allen_cahn] * stencil_phase.D
    rates += [1.0] * (stencil_phase.Q - stencil_phase.D - 1)

    lbm_config_phase = LBMConfig(stencil=stencil_phase, method=Method.MRT, compressible=True,
                                 delta_equilibrium=False,
                                 force=sp.symbols(f"F_:{stencil_phase.D}"), velocity_input=u,
                                 weighted=True, relaxation_rates=rates,
                                 output={'density': C_tmp}, kernel_type='stream_pull_collide')
    method_phase = create_lb_method(lbm_config=lbm_config_phase)

    lbm_config_hydro = LBMConfig(stencil=stencil_hydro, method=Method.MRT, compressible=False,
                                 weighted=True, relaxation_rate=omega,
                                 force=sp.symbols(f"F_:{stencil_hydro.D}"),
                                 output={'velocity': u}, kernel_type='collide_stream_push')
    method_hydro = create_lb_method(lbm_config=lbm_config_hydro)

    # create the kernels for the initialization of the g and h field
    h_updates = initializer_kernel_phase_field_lb(method_phase, C, u, h, parameters)
    h_updates = h_updates.new_with_substitutions(parameters.parameter_map())

    g_updates = initializer_kernel_hydro_lb(method_hydro, 1, u, g)
    g_updates = g_updates.new_with_substitutions(parameters.parameter_map())

    force_h = interface_tracking_force(C, stencil_phase, parameters)
    hydro_force = hydrodynamic_force(C, method_hydro, parameters, body_force)

    ####################
    # LBM UPDATE RULES #
    ####################

    lbm_optimisation = LBMOptimisation(symbolic_field=h, symbolic_temporary_field=h_tmp)
    phase_field_LB_step = create_lb_update_rule(lbm_config=lbm_config_phase,
                                                lbm_optimisation=lbm_optimisation)

    phase_field_LB_step = add_interface_tracking_force(phase_field_LB_step, force_h)
    phase_field_LB_step = phase_field_LB_step.new_with_substitutions(parameters.parameter_map())

    phase_field_LB_step = sympy_cse(phase_field_LB_step)
    # ---------------------------------------------------------------------------------------------------------
    force_Assignments = hydrodynamic_force_assignments(u, C, method_hydro, parameters, body_force)

    lbm_optimisation = LBMOptimisation(symbolic_field=g, symbolic_temporary_field=g_tmp)
    hydro_LB_step = create_lb_update_rule(lbm_config=lbm_config_hydro,
                                          lbm_optimisation=lbm_optimisation)

    hydro_LB_step = add_hydrodynamic_force(hydro_LB_step, force_Assignments, C, g, parameters)
    hydro_LB_step = hydro_LB_step.new_with_substitutions(parameters.parameter_map())

    hydro_LB_step = sympy_cse(hydro_LB_step)

    ###################
    # GENERATE SWEEPS #
    ###################

    # by default NT Stores are deactivated because they do not work in all cases
    # must be activated to achieve full potential for example on AVX512 CPUs
    cpu_vec = {'assume_inner_stride_one': True, 'nontemporal': False}

    vp = [('int32_t', 'cudaBlockSize0'),
          ('int32_t', 'cudaBlockSize1'),
          ('int32_t', 'cudaBlockSize2')]

    sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                        TypedSymbol("cudaBlockSize1", np.int32),
                        TypedSymbol("cudaBlockSize2", np.int32))
    sweep_params = {'block_size': sweep_block_size}

    stencil_typedefs = {'Stencil_phase_T': stencil_phase,
                        'Stencil_hydro_T': stencil_hydro}
    field_typedefs = {'PdfField_phase_T': h,
                      'PdfField_hydro_T': g,
                      'VelocityField_T': u,
                      'PhaseField_T': C}

    additional_code = f"""
    const char * StencilNamePhase = "{stencil_phase.name}";
    const char * StencilNameHydro = "{stencil_hydro.name}";
    """

    if not ctx.cuda:
        if not ctx.optimize_for_localhost:
            cpu_vec = {'instruction_set': None}

        generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates, target=Target.CPU)
        generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates, target=Target.CPU)

        generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                       field_swaps=[(h, h_tmp), (C, C_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec,
                       target=Target.CPU)

        generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                       field_swaps=[(g, g_tmp)],
                       inner_outer_split=True,
                       cpu_vectorize_info=cpu_vec,
                       target=Target.CPU)

        # communication
        generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                              streaming_pattern='pull', target=Target.CPU)

        generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                              streaming_pattern='push', target=Target.CPU)

        generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target=Target.CPU)

    if ctx.cuda:
        generate_sweep(ctx, 'initialize_phase_field_distributions',
                       h_updates, target=Target.GPU)
        generate_sweep(ctx, 'initialize_velocity_based_distributions',
                       g_updates, target=Target.GPU)

        generate_sweep(ctx, 'phase_field_LB_step', phase_field_LB_step,
                       field_swaps=[(h, h_tmp), (C, C_tmp)],
                       target=Target.GPU,
                       gpu_indexing_params=sweep_params,
                       varying_parameters=vp)

        generate_sweep(ctx, 'hydro_LB_step', hydro_LB_step,
                       field_swaps=[(g, g_tmp)],
                       target=Target.GPU,
                       gpu_indexing_params=sweep_params,
                       varying_parameters=vp)
        # communication
        generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                              streaming_pattern='pull', target=Target.GPU)

        generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                              streaming_pattern='push', target=Target.GPU)

        generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target=Target.GPU)

        # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'GenDefines', stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code=additional_code)
