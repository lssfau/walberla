from lbmpy_walberla import generate_alternating_lbm_sweep, generate_boundary, generate_alternating_lbm_boundary
from lbmpy_walberla.additional_data_handler import OutflowAdditionalDataHandler
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_ast
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow
from lbmpy.advanced_streaming import Timestep

from pystencils import Field
from lbmpy.stencils import get_stencil

#   Common Setup

stencil = get_stencil('D3Q27')
dim = len(stencil[0])
q = len(stencil)
target = 'cpu'
inplace_pattern = 'aa'
two_fields_pattern = 'pull'
namespace = 'lbmpy'

f_field = Field.create_generic('f', dim, index_shape=(q,), layout='fzyx')
f_field_tmp = Field.create_generic('f_tmp', dim, index_shape=(q,), layout='fzyx')
u_field = Field.create_generic('u', dim, index_shape=(dim,), layout='fzyx')

output = {
    'velocity': u_field
}

method_params = {
    'stencil': stencil,
    'method': 'srt',
    'relaxation_rate': 1.5,
    'output': output
}

optimization = {
    'target': target,
    'symbolic_field': f_field,
    'symbolic_temporary_field': f_field_tmp
}

collision_rule = create_lb_collision_rule(**method_params)
lb_method = collision_rule.method
noslip = NoSlip()
ubb = UBB((0.05,) + (0,) * (dim - 1))

outflow_normal = (1,) + (0,) * (dim - 1)
outflow_pull = ExtrapolationOutflow(outflow_normal, lb_method, streaming_pattern=two_fields_pattern)

outflow_inplace = ExtrapolationOutflow(outflow_normal, lb_method, streaming_pattern=inplace_pattern)

init_velocity = (0, ) * dim

init_kernel_pull = macroscopic_values_setter(lb_method, 1, init_velocity, f_field, streaming_pattern=two_fields_pattern)
init_kernel_inplace = macroscopic_values_setter(
    lb_method, 1, init_velocity, f_field, streaming_pattern=inplace_pattern, previous_timestep=Timestep.ODD)

stencil_typedefs = {'Stencil_T': stencil}
field_typedefs = {'PdfField_T': f_field, 'VelocityField_T': u_field}

with CodeGeneration() as ctx:
    #   Pull-Pattern classes
    ast_pull = create_lb_ast(collision_rule=collision_rule,
                             streaming_pattern=two_fields_pattern, optimization=optimization)
    generate_sweep(ctx, 'PullSweep', ast_pull, field_swaps=[(f_field, f_field_tmp)], namespace=namespace)

    generate_boundary(ctx, 'PullNoSlip', noslip, lb_method,
                      streaming_pattern=two_fields_pattern, target=target, namespace=namespace)
    generate_boundary(ctx, 'PullUBB', ubb, lb_method, streaming_pattern=two_fields_pattern,
                      target=target, namespace=namespace)
    generate_boundary(ctx, 'PullOutflow', outflow_pull, lb_method,
                      streaming_pattern=two_fields_pattern, target=target, namespace=namespace)

    generate_sweep(ctx, 'PullInit', init_kernel_pull, target=target, namespace=namespace)

    #   Inplace Pattern classes
    generate_alternating_lbm_sweep(ctx, 'InPlaceSweep', collision_rule, inplace_pattern,
                                   optimization=optimization, namespace=namespace)

    generate_alternating_lbm_boundary(ctx, 'InPlaceNoSlip', noslip, lb_method, streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)
    generate_alternating_lbm_boundary(ctx, 'InPlaceUBB', ubb, lb_method, streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)
    generate_alternating_lbm_boundary(ctx, 'InPlaceOutflow', outflow_inplace, lb_method, streaming_pattern=inplace_pattern,
                                      after_collision=True, target=target, namespace=namespace)

    generate_sweep(ctx, 'InPlaceInit', init_kernel_inplace, target=target, namespace=namespace)

    generate_info_header(ctx, "InplaceStreamingCodegen.h",
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
