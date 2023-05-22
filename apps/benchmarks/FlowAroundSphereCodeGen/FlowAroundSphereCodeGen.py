from pystencils import Target
from pystencils.field import fields
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil

from lbmpy.advanced_streaming.utility import get_timesteps
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow

from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header

from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler, OutflowAdditionalDataHandler
from lbmpy_walberla import generate_lb_pack_info
from lbmpy_walberla import generate_alternating_lbm_sweep, generate_alternating_lbm_boundary

import sympy as sp

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)
    q = stencil.Q
    dim = stencil.D
    streaming_pattern = 'esotwist'
    timesteps = get_timesteps(streaming_pattern)

    pdfs, velocity_field, density_field = fields(f"pdfs({q}), velocity({dim}), density(1) : {data_type}[{dim}D]",
                                                 layout='fzyx')
    omega = sp.Symbol("omega")
    u_max = sp.Symbol("u_max")

    output = {
        'density': density_field,
        'velocity': velocity_field
    }

    lbm_config = LBMConfig(stencil=stencil, method=Method.CUMULANT, compressible=True,
                           relaxation_rate=omega, galilean_correction=True,
                           field_name='pdfs', streaming_pattern=streaming_pattern, output=output)

    lbm_optimisation = LBMOptimisation(symbolic_field=pdfs, cse_global=False, cse_pdfs=False)

    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_optimisation)
    lb_method = collision_rule.method

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs, density=1.0,
                                                   streaming_pattern=streaming_pattern,
                                                   previous_timestep=timesteps[0])
    setter_assignments = setter_assignments.new_without_unused_subexpressions()

    # opt = {'instruction_set': 'sse', 'assume_aligned': True, 'nontemporal': False, 'assume_inner_stride_one': True}

    stencil_typedefs = {'Stencil_T': stencil}
    field_typedefs = {'PdfField_T': pdfs,
                      'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    if ctx.cuda:
        target = Target.GPU
    else:
        target = Target.CPU

    # sweeps
    generate_alternating_lbm_sweep(ctx, 'FlowAroundSphereCodeGen_LbSweep',
                                   collision_rule, lbm_config=lbm_config, lbm_optimisation=lbm_optimisation,
                                   target=target)
    generate_sweep(ctx, 'FlowAroundSphereCodeGen_MacroSetter', setter_assignments, target=target)

    # boundaries
    ubb = UBB(lambda *args: None, dim=dim, data_type=data_type)
    ubb_data_handler = UBBAdditionalDataHandler(stencil, ubb)
    outflow = ExtrapolationOutflow(stencil[4], lb_method, streaming_pattern=streaming_pattern, data_type=data_type)
    outflow_data_handler = OutflowAdditionalDataHandler(stencil, outflow, target=target)

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_UBB', ubb, lb_method,
                                      target=target, streaming_pattern=streaming_pattern,
                                      additional_data_handler=ubb_data_handler,
                                      layout='fzyx')

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_NoSlip', NoSlip(), lb_method,
                                      target=target, streaming_pattern=streaming_pattern,
                                      layout='fzyx')

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_Outflow', outflow, lb_method,
                                      target=target, streaming_pattern=streaming_pattern,
                                      additional_data_handler=outflow_data_handler,
                                      layout='fzyx')

    # communication
    generate_lb_pack_info(ctx, 'FlowAroundSphereCodeGen_PackInfo', stencil, pdfs,
                          streaming_pattern=streaming_pattern, always_generate_separate_classes=True, target=target)

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'FlowAroundSphereCodeGen_InfoHeader',
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
