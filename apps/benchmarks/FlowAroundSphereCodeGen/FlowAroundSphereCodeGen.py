from pystencils.field import fields

from lbmpy.advanced_streaming.utility import get_timesteps, Timestep
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.stencils import get_stencil
from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow

from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header

from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler, OutflowAdditionalDataHandler
from lbmpy_walberla import generate_boundary, generate_lb_pack_info
from lbmpy_walberla import generate_alternating_lbm_sweep, generate_alternating_lbm_boundary

import sympy as sp

stencil = get_stencil("D3Q27")
q = len(stencil)
dim = len(stencil[0])
streaming_pattern = 'esotwist'
timesteps = get_timesteps(streaming_pattern)

pdfs, velocity_field, density_field = fields(f"pdfs({q}), velocity({dim}), density(1) : double[{dim}D]", layout='fzyx')
omega = sp.Symbol("omega")
u_max = sp.Symbol("u_max")

output = {
    'density': density_field,
    'velocity': velocity_field
}

opt = {'symbolic_field': pdfs,
       'cse_global': False,
       'cse_pdfs': False}

method_params = {'method': 'cumulant',
                 'stencil': stencil,
                 'relaxation_rate': omega,
                 'galilean_correction': True,
                 'field_name': 'pdfs',
                 'streaming_pattern': streaming_pattern,
                 'output': output,
                 'optimization': opt}

collision_rule = create_lb_collision_rule(**method_params)
lb_method = collision_rule.method

# getter & setter
setter_assignments = macroscopic_values_setter(lb_method, velocity=velocity_field.center_vector,
                                               pdfs=pdfs, density=1,
                                               streaming_pattern=streaming_pattern,
                                               previous_timestep=timesteps[0])

# opt = {'instruction_set': 'sse', 'assume_aligned': True, 'nontemporal': False, 'assume_inner_stride_one': True}

stencil_typedefs = {'Stencil_T': stencil}
field_typedefs = {'PdfField_T': pdfs,
                  'VelocityField_T': velocity_field,
                  'ScalarField_T': density_field}

with CodeGeneration() as ctx:
    if ctx.cuda:
        target = 'gpu'
    else:
        target = 'cpu'

    opt['target'] = target

    # sweeps
    generate_alternating_lbm_sweep(ctx, 'FlowAroundSphereCodeGen_LbSweep',
                                   collision_rule, streaming_pattern, optimization=opt)
    generate_sweep(ctx, 'FlowAroundSphereCodeGen_MacroSetter', setter_assignments, target=target)

    # boundaries
    ubb = UBB(lambda *args: None, dim=dim)
    ubb_data_handler = UBBAdditionalDataHandler(stencil, ubb)
    outflow = ExtrapolationOutflow(stencil[4], lb_method, streaming_pattern=streaming_pattern)
    outflow_data_handler = OutflowAdditionalDataHandler(stencil, outflow, target=target)

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_UBB', ubb, lb_method,
                                      target=target, streaming_pattern=streaming_pattern,
                                      additional_data_handler=ubb_data_handler)

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_NoSlip', NoSlip(), lb_method,
                                      target=target, streaming_pattern=streaming_pattern)

    generate_alternating_lbm_boundary(ctx, 'FlowAroundSphereCodeGen_Outflow', outflow, lb_method,
                                      target=target, streaming_pattern=streaming_pattern,
                                      additional_data_handler=outflow_data_handler)

    # communication
    generate_lb_pack_info(ctx, 'FlowAroundSphereCodeGen_PackInfo', stencil, pdfs,
                          streaming_pattern=streaming_pattern, always_generate_separate_classes=True, target=target)

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'FlowAroundSphereCodeGen_InfoHeader',
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
