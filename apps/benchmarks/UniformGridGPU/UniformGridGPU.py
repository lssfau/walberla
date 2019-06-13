import sympy as sp
import numpy as np
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB
from pystencils_walberla import generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model, generate_boundary
from pystencils_walberla import CodeGeneration, generate_sweep
from pystencils.data_types import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions

omega = sp.symbols("omega")
# sweep_block_size = (128, 1, 1)
sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                    TypedSymbol("cudaBlockSize1", np.int32),
                    1)

sweep_params = {'block_size': sweep_block_size}

options = {
    'srt': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'relaxation_rate': omega,
        'compressible': False,
    },
    'trt': {
        'method': 'trt',
        'stencil': 'D3Q19',
        'relaxation_rate': omega,
    },
    'mrt': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [0, omega, 1.3, 1.4, omega, 1.2, 1.1],
    },
    'entropic': {
        'method': 'mrt3',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rates': [omega, omega, sp.Symbol("omega_free")],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'smagorinsky': True,
        'relaxation_rate': omega,
    }
}

with CodeGeneration() as ctx:

    common_options = {
        'field_name': 'pdfs',
        'temporary_field_name': 'pdfs_tmp',
        #'kernel_type': 'collide_stream_push',
        'optimization': {'cse_global': True,
                         'cse_pdfs': False}
    }
    options = options[ctx.config]
    options.update(common_options)

    vp = [
        ('int32_t', 'cudaBlockSize0'),
        ('int32_t', 'cudaBlockSize1')
    ]
    lb_method = create_lb_method(**options)
    update_rule = create_lb_update_rule(lb_method=lb_method, **options)

    update_rule = insert_fast_divisions(update_rule)
    update_rule = insert_fast_sqrts(update_rule)

    # CPU lattice model - required for macroscopic value computation, VTK output etc.
    options_without_opt = options.copy()
    del options_without_opt['optimization']
    generate_lattice_model(ctx, 'UniformGridGPU_LatticeModel', lb_method, update_rule_params=options_without_opt)

    # gpu LB sweep & boundaries
    generate_sweep(ctx, 'UniformGridGPU_LbKernel', update_rule,
                   field_swaps=[('pdfs', 'pdfs_tmp')],
                   inner_outer_split=True, target='gpu', gpu_indexing_params=sweep_params,
                   varying_parameters=vp)

    generate_boundary(ctx, 'UniformGridGPU_NoSlip', NoSlip(), lb_method, target='gpu')
    generate_boundary(ctx, 'UniformGridGPU_UBB', UBB([0.05, 0, 0]), lb_method, target='gpu')

    # communication
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_PackInfo', update_rule, target='gpu')
