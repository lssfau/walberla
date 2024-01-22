import numpy as np
import sympy as sp

from pystencils import fields, Target
from pystencils.astnodes import Conditional, Block, SympyAssignment
from pystencils.node_collection import NodeCollection
from pystencils.typing import TypedSymbol
from pystencils.simp.simplifications import sympy_cse

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, create_lb_update_rule
from lbmpy.creationfunctions import create_lb_method
from lbmpy.boundaries import NoSlip

from lbmpy.phasefield_allen_cahn.contact_angle import ContactAngle
from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, hydrodynamic_force, add_interface_tracking_force, \
    add_hydrodynamic_force, hydrodynamic_force_assignments
from lbmpy.phasefield_allen_cahn.parameter_calculation import AllenCahnParameters

import pystencils_walberla
from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_for_field, generate_info_header
from lbmpy_walberla import generate_boundary, generate_lb_pack_info

with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"

    contact_angle_in_degrees = 22

    stencil_phase = LBStencil(Stencil.D3Q27)
    stencil_hydro = LBStencil(Stencil.D3Q27)
    assert (stencil_phase.D == stencil_hydro.D)

    # In the codegeneration skript we only need the symbolic parameters
    parameters = AllenCahnParameters(0, 0, 0, 0, 0)

    ########################
    # FIELDS #
    ########################

    # velocity field
    u = fields(f"vel_field({stencil_hydro.D}): {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    # phase-field
    C = fields(f"phase_field: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    C_tmp = fields(f"phase_field_tmp: {field_type}[{stencil_hydro.D}D]", layout='fzyx')

    flag = fields(f"flag_field: uint8[{stencil_hydro.D}D]", layout='fzyx')
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
    g_updates = initializer_kernel_hydro_lb(method_hydro, 1, u, g)

    force_h = interface_tracking_force(C, stencil_phase, parameters)
    hydro_force = hydrodynamic_force(C, method_hydro, parameters, body_force)

    ####################
    # LBM UPDATE RULES #
    ####################

    lbm_optimisation = LBMOptimisation(symbolic_field=h, symbolic_temporary_field=h_tmp)
    allen_cahn_update_rule = create_lb_update_rule(lbm_config=lbm_config_phase,
                                                   lbm_optimisation=lbm_optimisation)

    allen_cahn_update_rule = add_interface_tracking_force(allen_cahn_update_rule, force_h)

    allen_cahn_update_rule = sympy_cse(allen_cahn_update_rule)

    allen_cahn_update_rule = NodeCollection([Conditional(sp.Eq(flag.center(), 2),
                                             Block(allen_cahn_update_rule.all_assignments),
                                             Block([SympyAssignment(C_tmp.center, C.center)]))])
    # ---------------------------------------------------------------------------------------------------------
    force_Assignments = hydrodynamic_force_assignments(u, C, method_hydro, parameters, body_force)

    lbm_optimisation = LBMOptimisation(symbolic_field=g, symbolic_temporary_field=g_tmp)
    hydro_lb_update_rule = create_lb_update_rule(lbm_config=lbm_config_hydro,
                                                 lbm_optimisation=lbm_optimisation)

    hydro_lb_update_rule = add_hydrodynamic_force(hydro_lb_update_rule, force_Assignments, C, g, parameters)

    hydro_lb_update_rule = NodeCollection([Conditional(sp.Eq(flag.center(), 2),
                                           Block(hydro_lb_update_rule.all_assignments))])

    contact_angle = ContactAngle(contact_angle_in_degrees, parameters.symbolic_interface_thickness,
                                 data_type=field_type)

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

    generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates, target=Target.GPU)
    generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates, target=Target.GPU)

    generate_sweep(ctx, 'phase_field_LB_step', allen_cahn_update_rule,
                   field_swaps=[(h, h_tmp), (C, C_tmp)],
                   target=Target.GPU,
                   gpu_indexing_params=sweep_params,
                   varying_parameters=vp)
    generate_boundary(ctx, 'phase_field_LB_NoSlip', NoSlip(), method_phase, target=Target.GPU, streaming_pattern='pull')

    generate_sweep(ctx, 'hydro_LB_step', hydro_lb_update_rule,
                   field_swaps=[(g, g_tmp)],
                   target=Target.GPU,
                   gpu_indexing_params=sweep_params,
                   varying_parameters=vp)
    generate_boundary(ctx, 'hydro_LB_NoSlip', NoSlip(), method_hydro, target=Target.GPU, streaming_pattern='push')

    # communication

    generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                          streaming_pattern='pull', target=Target.GPU)

    generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                          streaming_pattern='push', target=Target.GPU)

    generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target=Target.GPU)

    pystencils_walberla.boundary.generate_boundary(ctx, 'ContactAngle', contact_angle,
                                                   C.name, stencil_hydro, index_shape=[], target=Target.GPU)

    generate_info_header(ctx, 'GenDefines', stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code=additional_code)
