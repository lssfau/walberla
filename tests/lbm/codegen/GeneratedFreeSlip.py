from dataclasses import replace

from pystencils.field import fields
from lbmpy.advanced_streaming.utility import get_timesteps, is_inplace
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Stencil, Method
from lbmpy.creationfunctions import create_lb_method, create_lb_collision_rule
from lbmpy.boundaries import NoSlip, FreeSlip
from lbmpy_walberla.additional_data_handler import FreeSlipAdditionalDataHandler
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header
from lbmpy_walberla import generate_boundary, generate_lb_pack_info, generate_alternating_lbm_boundary, generate_alternating_lbm_sweep

import sympy as sp

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)

    pdfs, pdfs_tmp = fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]",
                            layout='fzyx')
    velocity_field, density_field = fields(f"velocity({stencil.D}), density(1) : {data_type}[{stencil.D}D]",
                                           layout='fzyx')
    omega = sp.Symbol("omega")

    output = {
        'density': density_field,
        'velocity': velocity_field
    }

    streaming_pattern = 'esotwist'
    timesteps = get_timesteps(streaming_pattern)

    lbm_config = LBMConfig(method=Method.SRT, stencil=stencil, relaxation_rate=omega, force=(1e-5, 0, 0),
                           output=output, streaming_pattern=streaming_pattern)

    lbm_opt = LBMOptimisation(symbolic_field=pdfs,
                              cse_global=False, cse_pdfs=False)

    method = create_lb_method(lbm_config=lbm_config)

    # getter & setter
    setter_assignments = macroscopic_values_setter(method, density=density_field.center,
                                                   velocity=velocity_field.center_vector,
                                                   pdfs=pdfs, streaming_pattern=streaming_pattern,
                                                   previous_timestep=timesteps[0])

    if not is_inplace(streaming_pattern):
        lbm_opt = replace(lbm_opt, symbolic_temporary_field=pdfs_tmp)
        field_swaps = [(pdfs, pdfs_tmp)]
    else:
        field_swaps = []

    collision_rule = create_lb_collision_rule(lb_method=method, lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    stencil_typedefs = {'Stencil_T': stencil}
    field_typedefs = {'PdfField_T': pdfs,
                      'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    # sweeps
    generate_alternating_lbm_sweep(ctx, 'GeneratedFreeSlip_Sweep', collision_rule, lbm_config=lbm_config,
                                   lbm_optimisation=lbm_opt, field_swaps=field_swaps)
    generate_sweep(ctx, 'GeneratedFreeSlip_MacroSetter', setter_assignments)

    # boundaries
    generate_alternating_lbm_boundary(ctx, 'GeneratedFreeSlip_NoSlip', NoSlip(), method,
                                      streaming_pattern=streaming_pattern)

    free_slip = FreeSlip(stencil=stencil)
    free_slip_data_handler = FreeSlipAdditionalDataHandler(stencil, free_slip)

    generate_alternating_lbm_boundary(ctx, 'GeneratedFreeSlip_FreeSlip', free_slip, method,
                                      additional_data_handler=free_slip_data_handler,
                                      streaming_pattern=streaming_pattern)

    # communication
    generate_lb_pack_info(ctx, 'GeneratedFreeSlip_PackInfo', stencil, pdfs, streaming_pattern=streaming_pattern,
                          always_generate_separate_classes=True)

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, "GeneratedFreeSlip.h",
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
