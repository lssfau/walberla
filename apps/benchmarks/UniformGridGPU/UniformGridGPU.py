import sympy as sp
import numpy as np
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.fieldaccess import StreamPullTwoFieldsAccessor, StreamPushTwoFieldsAccessor
from pystencils_walberla import generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model, generate_boundary
from pystencils_walberla import CodeGeneration, generate_sweep
from pystencils.data_types import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter, macroscopic_values_setter

omega = sp.symbols("omega")
compile_time_block_size = False

if compile_time_block_size:
    sweep_block_size = (128, 1, 1)
else:
    sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                        TypedSymbol("cudaBlockSize1", np.int32),
                        1)

sweep_params = {'block_size': sweep_block_size}

options_dict = {
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
    accessor = StreamPullTwoFieldsAccessor()
    #accessor = StreamPushTwoFieldsAccessor()
    assert not accessor.is_inplace, "This app does not work for inplace accessors"

    common_options = {
        'field_name': 'pdfs',
        'temporary_field_name': 'pdfs_tmp',
        'kernel_type': accessor,
        'optimization': {'cse_global': True,
                         'cse_pdfs': False}
    }
    options = options_dict.get(ctx.config, options_dict['srt'])
    options.update(common_options)

    stencil_str = options['stencil']
    q = int(stencil_str[stencil_str.find('Q')+1:])
    pdfs, velocity_field = ps.fields("pdfs({q}), velocity(3) : double[3D]".format(q=q), layout='fzyx')
    options['optimization']['symbolic_field'] = pdfs

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

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=1)
    getter_assignments = macroscopic_values_getter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector,  density=None)
    generate_sweep(ctx, 'UniformGridGPU_MacroSetter', setter_assignments)
    generate_sweep(ctx, 'UniformGridGPU_MacroGetter', getter_assignments)

    # communication
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_PackInfo', update_rule, target='gpu')
