from dataclasses import replace

from lbmpy_walberla import generate_alternating_lbm_sweep, generate_boundary, generate_alternating_lbm_boundary
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header

from pystencils import Target, CreateKernelConfig
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_ast
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow
from lbmpy.advanced_streaming import Timestep

from pystencils import Field

#   Common Setup

stencil = LBStencil(Stencil.D3Q27)
target = Target.CPU
inplace_pattern = 'aa'
two_fields_pattern = 'pull'
namespace = 'lbmpy'

f_field = Field.create_generic('f', stencil.D, index_shape=(stencil.Q,), layout='fzyx')
f_field_tmp = Field.create_generic('f_tmp', stencil.D, index_shape=(stencil.Q,), layout='fzyx')
u_field = Field.create_generic('u', stencil.D, index_shape=(stencil.D,), layout='fzyx')

output = {
    'velocity': u_field
}

lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=1.5, output=output)
lbm_opt = LBMOptimisation(symbolic_field=f_field,
                          symbolic_temporary_field=f_field_tmp)

config = CreateKernelConfig(target=target)

collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt, config=config)
lb_method = collision_rule.method
noslip = NoSlip()
ubb = UBB((0.05,) + (0,) * (stencil.D - 1))

outflow_normal = (1,) + (0,) * (stencil.D - 1)
outflow_pull = ExtrapolationOutflow(outflow_normal, lb_method, streaming_pattern=two_fields_pattern)

outflow_inplace = ExtrapolationOutflow(outflow_normal, lb_method, streaming_pattern=inplace_pattern)

init_velocity = (0,) * stencil.D

init_kernel_pull = macroscopic_values_setter(lb_method, 1, init_velocity, f_field, streaming_pattern=two_fields_pattern)
init_kernel_inplace = macroscopic_values_setter(
    lb_method, 1, init_velocity, f_field, streaming_pattern=inplace_pattern, previous_timestep=Timestep.ODD)

stencil_typedefs = {'Stencil_T': stencil}
field_typedefs = {'PdfField_T': f_field, 'VelocityField_T': u_field}

with CodeGeneration() as ctx:
    #   Pull-Pattern classes
    ast_pull = create_lb_ast(collision_rule=collision_rule,
                             streaming_pattern=two_fields_pattern, lbm_optimisation=lbm_opt)
    generate_sweep(ctx, 'PullSweep', ast_pull, field_swaps=[(f_field, f_field_tmp)], namespace=namespace)

    generate_boundary(ctx, 'PullNoSlip', noslip, lb_method,
                      streaming_pattern=two_fields_pattern, target=target, namespace=namespace)
    generate_boundary(ctx, 'PullUBB', ubb, lb_method, streaming_pattern=two_fields_pattern,
                      target=target, namespace=namespace)
    generate_boundary(ctx, 'PullOutflow', outflow_pull, lb_method,
                      streaming_pattern=two_fields_pattern, target=target, namespace=namespace)

    generate_sweep(ctx, 'PullInit', init_kernel_pull, target=target, namespace=namespace)

    #   Inplace Pattern classes
    inplace_lbm_config = replace(lbm_config, streaming_pattern=inplace_pattern)
    generate_alternating_lbm_sweep(ctx, 'InPlaceSweep', collision_rule,
                                   lbm_config=inplace_lbm_config, namespace=namespace)

    generate_alternating_lbm_boundary(ctx, 'InPlaceNoSlip', noslip, lb_method, streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)
    generate_alternating_lbm_boundary(ctx, 'InPlaceUBB', ubb, lb_method, streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)
    generate_alternating_lbm_boundary(ctx, 'InPlaceOutflow', outflow_inplace, lb_method,
                                      streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)

    generate_sweep(ctx, 'InPlaceInit', init_kernel_inplace, target=target, namespace=namespace)

    generate_info_header(ctx, "InplaceStreamingCodegen.h",
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
