from dataclasses import replace

from pystencils.field import fields
from lbmpy.advanced_streaming.utility import get_timesteps, is_inplace
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Stencil, Method, lattice_viscosity_from_relaxation_rate
from lbmpy.creationfunctions import create_lb_method, create_lb_collision_rule
from lbmpy.boundaries import NoSlip, FreeSlip, WallFunctionBounce, SpaldingsLaw
from lbmpy_walberla.additional_data_handler import FreeSlipAdditionalDataHandler
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header
from lbmpy_walberla import generate_boundary, generate_lb_pack_info, generate_alternating_lbm_boundary, generate_alternating_lbm_sweep

import sympy as sp

from lbmpy_walberla.additional_data_handler import WFBAdditionalDataHandler

with CodeGeneration() as ctx:

    config_tokens = ctx.config.split('_')

    assert len(config_tokens) == 3
    stencil_str = config_tokens[0]
    wfb_weights_str = config_tokens[1]
    velocity_str = config_tokens[2]

    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)

    pdfs, pdfs_tmp = fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]",
                            layout='fzyx')
    velocity_field = fields(f"velocity({stencil.D}): {data_type}[{stencil.D}D]", layout='fzyx')
    omega = sp.Symbol("omega")

    output = {
        'velocity': velocity_field
    }

    streaming_pattern = 'aa'
    timesteps = get_timesteps(streaming_pattern)

    lbm_config = LBMConfig(method=Method.SRT, stencil=stencil, relaxation_rate=omega, force=(1e-5, 0, 0),
                           output=output, streaming_pattern=streaming_pattern)

    lbm_opt = LBMOptimisation(symbolic_field=pdfs,
                              cse_global=False, cse_pdfs=False)

    lbm_method = create_lb_method(lbm_config=lbm_config)

    if not is_inplace(streaming_pattern):
        lbm_opt = replace(lbm_opt, symbolic_temporary_field=pdfs_tmp)

    collision_rule = create_lb_collision_rule(lb_method=lbm_method, lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    stencil_typedefs = {'Stencil_T': stencil}
    field_typedefs = {'PdfField_T': pdfs,
                      'VelocityField_T': velocity_field}

    #   Boundary conditions
    nu = lattice_viscosity_from_relaxation_rate(omega)
    u_tau_target = sp.Symbol("target_u_tau")

    noslip = NoSlip()

    wall_function_model = SpaldingsLaw(viscosity=nu, kappa=0.41, b=5.5, newton_steps=5)
    wfb_common_params = {
        'lb_method': lbm_method,
        'velocity_field': velocity_field,
        'wall_function_model': wall_function_model,
        'use_maronga_correction': False,
        'sampling_shift': 0,
        'filter_width': sp.Symbol("filter_width"),
        'data_type': data_type,
        'target_friction_velocity': u_tau_target
    }

    if wfb_weights_str == 'lattice':
        wfb_common_params['weight_method'] = WallFunctionBounce.WeightMethod.LATTICE_WEIGHT
    elif wfb_weights_str == 'geometric':
        wfb_common_params['weight_method'] = WallFunctionBounce.WeightMethod.GEOMETRIC_WEIGHT
    elif wfb_weights_str == 'none':
        pass
    else:
        raise ValueError("Invalid weight method")

    if velocity_str == 'mean':
        wfb_common_params['reference_velocity'] = WallFunctionBounce.ReferenceVelocity.MEAN_VELOCITY
    elif velocity_str == 'filtered':
        wfb_common_params['reference_velocity'] = WallFunctionBounce.ReferenceVelocity.FILTERED_VELOCITY
    else:
        raise ValueError("Invalid reference velocity")

    normal_direction = (0, 1, 0)
    wfb = WallFunctionBounce(**wfb_common_params, normal_direction=normal_direction)

    additional_data_handler = WFBAdditionalDataHandler(stencil, wfb, velocity_field)
    generate_boundary(ctx, "WFB", wfb, lbm_method, additional_data_handler=additional_data_handler)

    config_definitions = []
    if wfb_common_params['use_maronga_correction']:
        config_definitions.append("#define RUN_WITH_MARONGA_CORRECTION")
    if velocity_str == 'mean':
        config_definitions.append("#define MEAN_VELOCITY")
    elif velocity_str == 'filtered':
        config_definitions.append("#define FILTERED_VELOCITY")

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, "WFBHeader.h",
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code="\n".join(config_definitions))
