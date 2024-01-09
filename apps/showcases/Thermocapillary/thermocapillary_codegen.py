import logging
import numpy as np
import sympy as sp

from pystencils import Assignment, fields, Target
from pystencils.astnodes import LoopOverCoordinate
from pystencils.typing import TypedSymbol
from pystencils.simp import (add_subexpressions_for_field_reads, sympy_cse, insert_aliases, insert_constants,
                             insert_symbol_times_minus_one, insert_squares,
                             insert_constant_multiples, insert_constant_additions)

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, create_lb_update_rule
from lbmpy.creationfunctions import create_lb_method
from lbmpy.boundaries import DiffusionDirichlet, NeumannByCopy, NoSlip, UBB
from lbmpy.macroscopic_value_kernels import pdf_initialization_assignments

from lbmpy.phasefield_allen_cahn.contact_angle import ContactAngle
from lbmpy.phasefield_allen_cahn.kernel_equations import initializer_kernel_phase_field_lb, \
    initializer_kernel_hydro_lb, interface_tracking_force, hydrodynamic_force, add_interface_tracking_force, \
    add_hydrodynamic_force, hydrodynamic_force_assignments
from lbmpy.phasefield_allen_cahn.numerical_solver import get_runge_kutta_update_assignments
from lbmpy.phasefield_allen_cahn.parameter_calculation import ThermocapillaryParameters

import pystencils_walberla
from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_for_field, generate_info_header
from lbmpy_walberla.additional_data_handler import DiffusionDirichletAdditionalDataHandler
from lbmpy_walberla import generate_boundary, generate_lb_pack_info

# to disable type conversion warnings
logging.disable(logging.WARNING)

one = sp.Number(1)
two = sp.Number(2)
three = sp.Number(3)
half = sp.Rational(1, 2)


with CodeGeneration() as ctx:
    field_type = "float64"
    field_type_pdfs = "float64"

    subs = {sp.Symbol(f"xi_{ctr}"): TypedSymbol(f"xi_{ctr}", dtype=field_type) for ctr in range(0, 100)}

    contact_angle_in_degrees = sp.Symbol("alpha")
    generate_with_heat_source = False

    stencil_phase = LBStencil(Stencil.D3Q19)
    stencil_hydro = LBStencil(Stencil.D3Q19)
    stencil_thermal = LBStencil(Stencil.D3Q7)
    assert (stencil_phase.D == stencil_hydro.D == stencil_thermal.D)
    full_stencil = LBStencil(Stencil.D2Q9) if stencil_phase.D == 2 else LBStencil(Stencil.D3Q27)

    # In the codegeneration skript we only need the symbolic parameters
    parameters = ThermocapillaryParameters(0, 0, 0, 0, 0, 0, 0)

    ########################
    # FIELDS #
    ########################

    # velocity field
    u = fields(f"vel_field({stencil_hydro.D}): {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    # phase-field
    C = fields(f"phase_field: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    C_tmp = fields(f"phase_field_tmp: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    # temperature field
    temperature = fields(f"temperature_field: {field_type}[{stencil_hydro.D}D]", layout='fzyx')

    # phase-field distribution functions
    h = fields(f"lb_phase_field({stencil_phase.Q}): {field_type_pdfs}[{stencil_phase.D}D]", layout='fzyx')
    h_tmp = fields(f"lb_phase_field_tmp({stencil_phase.Q}): {field_type_pdfs}[{stencil_phase.D}D]", layout='fzyx')
    # hydrodynamic distribution functions
    g = fields(f"lb_velocity_field({stencil_hydro.Q}): {field_type_pdfs}[{stencil_hydro.D}D]", layout='fzyx')
    g_tmp = fields(f"lb_velocity_field_tmp({stencil_hydro.Q}): {field_type_pdfs}[{stencil_hydro.D}D]", layout='fzyx')
    # thermal PDFs
    f = fields(f"lb_thermal_field({stencil_thermal.Q}): {field_type_pdfs}[{stencil_hydro.D}D]", layout='fzyx')
    f_tmp = fields(f"lb_thermal_field_tmp({stencil_thermal.Q}): {field_type_pdfs}[{stencil_hydro.D}D]", layout='fzyx')

    RK1 = fields(f"Rk1: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    RK2 = fields(f"Rk2: {field_type}[{stencil_hydro.D}D]", layout='fzyx')
    RK3 = fields(f"Rk3: {field_type}[{stencil_hydro.D}D]", layout='fzyx')

    ########################################
    # RELAXATION RATES AND EXTERNAL FORCES #
    ########################################
    k_l = parameters.symbolic_heat_conductivity_light
    k_h = parameters.symbolic_heat_conductivity_heavy

    conductivity = k_l + C.center * (k_h - k_l)
    relaxation_rate_thermal = one/(half + (three * conductivity))

    # relaxation rate for interface tracking LB step
    relaxation_rate_allen_cahn = one / (half + (three * parameters.symbolic_mobility))
    # calculate the relaxation rate for hydrodynamic LB step
    rho_L = parameters.symbolic_density_light
    rho_H = parameters.symbolic_density_heavy
    density = rho_L + C.center * (rho_H - rho_L)
    # force acting on all phases of the model
    body_force = [0, 0, 0]

    omega = parameters.omega(C)

    ###############
    # LBM METHODS #
    ###############

    lbm_config_phase = LBMConfig(stencil=stencil_phase, method=Method.CENTRAL_MOMENT, compressible=True,
                                 zero_centered=False,
                                 force=sp.symbols(f"F_:{stencil_phase.D}"), velocity_input=u,
                                 weighted=True, relaxation_rates=[relaxation_rate_allen_cahn, ] * stencil_phase.Q,
                                 output={'density': C_tmp})
    method_phase = create_lb_method(lbm_config=lbm_config_phase)

    lbm_config_hydro = LBMConfig(stencil=stencil_hydro, method=Method.CENTRAL_MOMENT, compressible=False,
                                 weighted=True, relaxation_rates=[omega, ] * stencil_hydro.Q,
                                 force=sp.symbols(f"F_:{stencil_hydro.D}"),
                                 output={'velocity': u})
    method_hydro = create_lb_method(lbm_config=lbm_config_hydro)

    config_thermal = LBMConfig(stencil=stencil_thermal, method=Method.CENTRAL_MOMENT, compressible=True,
                               zero_centered=False,
                               relaxation_rates=[relaxation_rate_thermal, ] * stencil_thermal.Q,
                               velocity_input=u, output={'density': temperature})
    method_thermal = create_lb_method(lbm_config=config_thermal)

    # create the kernels for the initialization of the g and h field
    h_updates = initializer_kernel_phase_field_lb(method_phase, C, u, h, parameters)
    g_updates = initializer_kernel_hydro_lb(method_hydro, 1.0, u, g)
    f_updates = pdf_initialization_assignments(method_thermal, temperature.center,
                                               u.center_vector, f.center_vector)

    force_h = interface_tracking_force(C, stencil_phase, parameters)
    hydro_force = hydrodynamic_force(C, method_hydro, parameters, body_force, temperature_field=temperature)

    heat_terms = list()
    counters = [LoopOverCoordinate.get_loop_counter_symbol(i) for i in range(stencil_thermal.D)]
    if generate_with_heat_source:
        block_offset = [TypedSymbol(f"block_offset_{i}", counter.dtype, nonnegative=True)
                        for i, counter in enumerate(counters)]
        ctr = [offset + counter for offset, counter in zip(counters, block_offset)]

        heat_source_midpoint_one = [TypedSymbol(f"hm_one_{i}", counter.dtype, nonnegative=True)
                                    for i, counter in enumerate(counters)]

        heat_source_midpoint_two = [TypedSymbol(f"hm_two_{i}", counter.dtype, nonnegative=True)
                                    for i, counter in enumerate(counters)]

        ws = sp.Symbol("w_s")
        size_diffused_hotspot = sp.Symbol("d_s")

        maximum_heat_flux_one = sp.Symbol("qs_one")
        maximum_heat_flux_two = sp.Symbol("qs_two")

        nominator_one = sum([(c - hm)**2 for c, hm in zip(ctr, heat_source_midpoint_one)])
        nominator_two = sum([(c - hm)**2 for c, hm in zip(ctr, heat_source_midpoint_two)])

        term_one = maximum_heat_flux_one * sp.exp(-2 * nominator_one / (ws**2))
        term_two = maximum_heat_flux_two * sp.exp(-2 * nominator_two / (ws**2))

        cond_one = sp.Piecewise((term_one, nominator_one <= size_diffused_hotspot**2), (0.0, True))
        cond_two = sp.Piecewise((term_two, nominator_two <= size_diffused_hotspot**2), (0.0, True))
        heat_source = cond_one + cond_two

        weights = method_thermal.weights

        for i in range(stencil_thermal.Q):
            heat_terms.append(weights[i] * heat_source)

    ####################
    # LBM UPDATE RULES #
    ####################

    lbm_optimisation = LBMOptimisation(symbolic_field=h, symbolic_temporary_field=h_tmp)
    allen_cahn_update_rule = create_lb_update_rule(lbm_config=lbm_config_phase,
                                                   lbm_optimisation=lbm_optimisation)

    allen_cahn_update_rule = add_interface_tracking_force(allen_cahn_update_rule, force_h)
    allen_cahn_update_rule = add_subexpressions_for_field_reads(allen_cahn_update_rule)
    allen_cahn_update_rule = allen_cahn_update_rule.new_with_substitutions(substitutions=subs)
    # ---------------------------------------------------------------------------------------------------------
    force_Assignments = hydrodynamic_force_assignments(u, C, method_hydro, parameters, body_force, sub_iterations=1,
                                                       temperature_field=temperature)

    lbm_optimisation = LBMOptimisation(symbolic_field=g, symbolic_temporary_field=g_tmp)
    hydro_lb_update_rule = create_lb_update_rule(lbm_config=lbm_config_hydro,
                                                 lbm_optimisation=lbm_optimisation)

    hydro_lb_update_rule = add_hydrodynamic_force(hydro_lb_update_rule, force_Assignments, C, g, parameters)
    hydro_lb_update_rule = sympy_cse(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_aliases(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_constants(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_symbol_times_minus_one(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_squares(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_constant_multiples(hydro_lb_update_rule)
    hydro_lb_update_rule = insert_constant_additions(hydro_lb_update_rule)
    hydro_lb_update_rule = add_subexpressions_for_field_reads(hydro_lb_update_rule)
    hydro_lb_update_rule = hydro_lb_update_rule.new_with_substitutions(substitutions=subs)
    # ---------------------------------------------------------------------------------------------------------
    lbm_optimisation = LBMOptimisation(symbolic_field=f, symbolic_temporary_field=f_tmp)
    thermal_lb_update_rule = create_lb_update_rule(lbm_config=config_thermal,
                                                   lbm_optimisation=lbm_optimisation)
    thermal_lb_update_rule = add_subexpressions_for_field_reads(thermal_lb_update_rule)
    thermal_lb_update_rule = thermal_lb_update_rule.new_with_substitutions(substitutions=subs)

    main_assignments = thermal_lb_update_rule.main_assignments

    if generate_with_heat_source:
        for i in range(stencil_thermal.Q):
            main_assignments[i] = Assignment(main_assignments[i].lhs, main_assignments[i].rhs + heat_terms[i])

    init_rk2 = [Assignment(RK1.center, temperature.center)]

    init_rk4 = [Assignment(RK1.center, temperature.center),
                Assignment(RK2.center, temperature.center),
                Assignment(RK3.center, temperature.center)]

    rk2_assignments = get_runge_kutta_update_assignments(full_stencil, C, temperature, u, [RK1, ],
                                                         conduction_h=k_h,
                                                         conduction_l=k_l,
                                                         heat_capacity_h=1,
                                                         heat_capacity_l=1,
                                                         density=density,
                                                         stabiliser=1)

    rk4_assignments = get_runge_kutta_update_assignments(full_stencil, C, temperature, u, [RK1, RK2, RK3],
                                                         conduction_h=k_h,
                                                         conduction_l=k_l,
                                                         heat_capacity_h=1,
                                                         heat_capacity_l=1,
                                                         density=density,
                                                         stabiliser=1)

    contact_angle = ContactAngle(contact_angle_in_degrees, parameters.symbolic_interface_thickness,
                                 data_type=field_type)

    ###################
    # GENERATE SWEEPS #
    ###################

    if ctx.gpu:
        target = Target.GPU
        openmp = False
        cpu_vec = None
        vp = [('int32_t', 'cudaBlockSize0'),
              ('int32_t', 'cudaBlockSize1'),
              ('int32_t', 'cudaBlockSize2')]

        sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                            TypedSymbol("cudaBlockSize1", np.int32),
                            TypedSymbol("cudaBlockSize2", np.int32))
        sweep_params = {'block_size': sweep_block_size}
    else:
        if ctx.optimize_for_localhost:
            cpu_vec = {"instruction_set": None}
        else:
            cpu_vec = None

        openmp = True if ctx.openmp else False

        target = Target.CPU
        sweep_params = {}
        vp = []

    comm_stencil_phase = LBStencil(Stencil.D3Q27) if stencil_phase == LBStencil(Stencil.D3Q15) else stencil_phase
    comm_stencil_hydro = LBStencil(Stencil.D3Q27) if stencil_hydro == LBStencil(Stencil.D3Q15) else stencil_hydro
    comm_stencil_thermal = LBStencil(Stencil.D3Q27) if stencil_thermal == LBStencil(Stencil.D3Q15) else stencil_thermal

    stencil_typedefs = {'Stencil_phase_T': stencil_phase, 'CommStencil_phase_T': comm_stencil_phase,
                        'Stencil_hydro_T': stencil_hydro, 'CommStencil_hydro_T': comm_stencil_hydro,
                        'Stencil_thermal_T': stencil_thermal, 'CommStencil_thermal_T': comm_stencil_thermal,
                        'Full_Stencil_T': full_stencil}
    field_typedefs = {'PdfField_phase_T': h,
                      'PdfField_hydro_T': g,
                      'PdfField_thermal_T': f,
                      'VectorField_T': u,
                      'ScalarField_T': C}

    additional_code = f"""
const char * StencilNamePhase = "{stencil_phase.name}";
const char * StencilNameHydro = "{stencil_hydro.name}";
const char * StencilNameThermal = "{stencil_thermal.name}";
const char * CollisionSpacePhase = "{lbm_config_phase.method.name.lower()}";
const char * CollisionSpaceHydro = "{lbm_config_hydro.method.name.lower()}";
const char * CollisionSpaceThermal = "{config_thermal.method.name.lower()}";

const char * fieldType = "{field_type}";
const char * fieldTypePDFs = "{field_type_pdfs}";
#if defined(WALBERLA_BUILD_WITH_GPU_SUPPORT)
using GPUFieldPDFs = {"walberla::gpu::GPUField< float >" 
    if field_type_pdfs == "float32" else "walberla::gpu::GPUField< double >"} ;
using GPUField = walberla::gpu::GPUField< walberla::real_t >;
#endif

{"#define GENERATED_HEAT_SOURCE" if generate_with_heat_source else ""}

#ifdef GENERATED_HEAT_SOURCE
const bool GeneratedHeatSource = true;
#else
const bool GeneratedHeatSource = false;
#endif
"""

    generate_sweep(ctx, 'initialize_phase_field_distributions', h_updates,
                   target=target, cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'initialize_velocity_based_distributions', g_updates,
                   target=target, cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'initialize_thermal_distributions', f_updates,
                   target=target, cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)

    generate_sweep(ctx, 'phase_field_LB_step', allen_cahn_update_rule,
                   field_swaps=[(h, h_tmp)],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_boundary(ctx, 'phase_field_LB_NoSlip', NoSlip(), method_phase,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)

    generate_sweep(ctx, 'hydro_LB_step', hydro_lb_update_rule,
                   field_swaps=[(g, g_tmp)],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec, max_threads=256)

    velocity_wall = tuple([TypedSymbol('w_x', dtype=field_type)] + [0.0] * (stencil_hydro.D - 1))
    generate_boundary(ctx, 'hydro_LB_NoSlip', NoSlip(), method_hydro,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)
    generate_boundary(ctx, 'hydro_LB_UBB', UBB(velocity_wall, data_type=field_type), method_hydro,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)

    if generate_with_heat_source:
        generate_sweep(ctx, 'thermal_LB_step', thermal_lb_update_rule,
                       field_swaps=[(f, f_tmp)],
                       target=target,
                       varying_parameters=vp, gpu_indexing_params=sweep_params,
                       cpu_openmp=openmp, cpu_vectorize_info=cpu_vec,
                       block_offset=block_offset)
    else:
        generate_sweep(ctx, 'thermal_LB_step', thermal_lb_update_rule,
                       field_swaps=[(f, f_tmp)],
                       target=target,
                       varying_parameters=vp, gpu_indexing_params=sweep_params,
                       cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)

    generate_sweep(ctx, 'initialize_rk2', init_rk2, ghost_layers_to_include=1,
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk2_first_stage', rk2_assignments[0],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk2_second_stage', rk2_assignments[1],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)

    generate_sweep(ctx, 'initialize_rk4', init_rk4,
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk4_first_stage', rk4_assignments[0],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk4_second_stage', rk4_assignments[1],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk4_third_stage', rk4_assignments[2],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)
    generate_sweep(ctx, 'rk4_fourth_stage', rk4_assignments[3],
                   target=target,
                   varying_parameters=vp, gpu_indexing_params=sweep_params,
                   cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)

    dirichlet_bc_static = DiffusionDirichlet(TypedSymbol('T_c', dtype=field_type), u, data_type=field_type)
    generate_boundary(ctx, 'thermal_LB_DiffusionDirichlet_static', dirichlet_bc_static, method_thermal,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)

    dirichlet_bc_dynamic = DiffusionDirichlet(lambda *args: None, u, data_type=field_type)
    diffusion_data_handler = DiffusionDirichletAdditionalDataHandler(stencil_thermal, dirichlet_bc_dynamic)
    generate_boundary(ctx, 'thermal_LB_DiffusionDirichlet_dynamic', dirichlet_bc_dynamic, method_thermal,
                      additional_data_handler=diffusion_data_handler,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)

    generate_boundary(ctx, 'thermal_LB_NeumannByCopy', NeumannByCopy(), method_thermal,
                      target=target, streaming_pattern='pull', cpu_openmp=openmp, data_type=field_type_pdfs)

    # communication

    generate_lb_pack_info(ctx, 'PackInfo_phase_field_distributions', stencil_phase, h,
                          streaming_pattern='pull', target=target)

    generate_lb_pack_info(ctx, 'PackInfo_velocity_based_distributions', stencil_hydro, g,
                          streaming_pattern='pull', target=target)

    generate_lb_pack_info(ctx, 'PackInfo_thermal_field_distributions', stencil_thermal, f,
                          streaming_pattern='pull', target=target)

    generate_pack_info_for_field(ctx, 'PackInfo_phase_field', C, target=target)
    generate_pack_info_for_field(ctx, 'PackInfo_temperature_field', temperature, target=target)

    pystencils_walberla.boundary.generate_boundary(ctx, 'ContactAngle', contact_angle,
                                                   C.name, stencil_hydro, index_shape=[], target=target)

    generate_info_header(ctx, 'GenDefines', stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code=additional_code)

