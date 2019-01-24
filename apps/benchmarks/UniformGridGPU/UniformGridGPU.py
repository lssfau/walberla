import sympy as sp
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB
from pystencils_walberla import generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model, generate_boundary
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    # LB options
    options = {
        'method': 'srt',
        'stencil': 'D3Q19',
        'relaxation_rate': sp.Symbol("omega"),
        'field_name': 'pdfs',
        'compressible': False,
        'temporary_field_name': 'pdfs_tmp',
        'optimization': {'cse_global': True,
                         'cse_pdfs': True,
                         'gpu_indexing_params': {'block_size': (128, 1, 1)}}
    }
    lb_method = create_lb_method(**options)
    update_rule = create_lb_update_rule(lb_method=lb_method, **options)

    # CPU lattice model - required for macroscopic value computation, VTK output etc.
    generate_lattice_model(ctx, 'UniformGridGPU_LatticeModel', lb_method)

    # gpu LB sweep & boundaries
    generate_sweep(ctx, 'UniformGridGPU_LbKernel', update_rule, field_swaps=[('pdfs', 'pdfs_tmp')],
                   inner_outer_split=True, target='gpu')
    generate_boundary(ctx, 'UniformGridGPU_NoSlip', NoSlip(), lb_method, target='gpu')
    generate_boundary(ctx, 'UniformGridGPU_UBB', UBB([0.05, 0, 0]), lb_method, target='gpu')

    # communication
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_PackInfo', update_rule, target='gpu')
