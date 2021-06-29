import sympy as sp
import numpy as np
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule, create_lb_collision_rule
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.fieldaccess import StreamPullTwoFieldsAccessor
from pystencils_walberla import generate_pack_info_from_kernel
from lbmpy_walberla import generate_lattice_model, generate_boundary
from pystencils_walberla import CodeGeneration, generate_sweep
from pystencils.data_types import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter, macroscopic_values_setter

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
omega_fill = sp.symbols("omega_:10")
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
        'relaxation_rates': [omega, 1.3, 1.4, 1.2, 1.1, 1.15, 1.234, 1.4235],
    },
    'mrt_full': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [omega_fill[0], omega, omega_fill[1], omega_fill[2],
                             omega_fill[3], omega_fill[4], omega_fill[5]],
    },
    'entropic': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rates': [omega, omega, omega_free, omega_free, omega_free, omega_free],
        'entropic': True,
    },
    'entropic_kbc_n4': {
        'method': 'trt-kbc-n4',
        'stencil': 'D3Q27',
        'compressible': True,
        'relaxation_rates': [omega, omega_free],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'smagorinsky': True,
        'relaxation_rate': omega,
    },
    'cumulant': {
        'method': 'cumulant',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rate': omega,
    },
}

info_header = """
#include "stencil/D3Q{q}.h"\nusing Stencil_T = walberla::stencil::D3Q{q};
const char * infoStencil = "{stencil}";
const char * infoConfigName = "{configName}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    accessor = StreamPullTwoFieldsAccessor()
    # accessor = StreamPushTwoFieldsAccessor()
    assert not accessor.is_inplace, "This app does not work for inplace accessors"

    common_options = {
        'field_name': 'pdfs',
        'temporary_field_name': 'pdfs_tmp',
        'kernel_type': accessor,
        'optimization': {'cse_global': True,
                         'cse_pdfs': False}
    }
    config_name = ctx.config
    noopt = False
    d3q27 = False
    if config_name.endswith("_noopt"):
        noopt = True
        config_name = config_name[:-len("_noopt")]
    if config_name.endswith("_d3q27"):
        d3q27 = True
        config_name = config_name[:-len("_d3q27")]

    options = options_dict[config_name]
    options.update(common_options)
    options = options.copy()

    if noopt:
        options['optimization']['cse_global'] = False
        options['optimization']['cse_pdfs'] = False
    if d3q27:
        options['stencil'] = 'D3Q27'

    stencil_str = options['stencil']
    q = int(stencil_str[stencil_str.find('Q') + 1:])
    pdfs, velocity_field = ps.fields("pdfs({q}), velocity(3) : double[3D]".format(q=q), layout='fzyx')
    options['optimization']['symbolic_field'] = pdfs

    vp = [
        ('double', 'omega_0'),
        ('double', 'omega_1'),
        ('double', 'omega_2'),
        ('double', 'omega_3'),
        ('double', 'omega_4'),
        ('double', 'omega_5'),
        ('double', 'omega_6'),
        ('int32_t', 'cudaBlockSize0'),
        ('int32_t', 'cudaBlockSize1'),
    ]
    lb_method = create_lb_method(**options)
    update_rule = create_lb_update_rule(lb_method=lb_method, **options)

    if not noopt:
        update_rule = insert_fast_divisions(update_rule)
        update_rule = insert_fast_sqrts(update_rule)

    # CPU lattice model - required for macroscopic value computation, VTK output etc.
    options_without_opt = options.copy()
    del options_without_opt['optimization']
    generate_lattice_model(ctx, 'UniformGridGPU_LatticeModel', create_lb_collision_rule(lb_method=lb_method,
                                                                                        **options_without_opt))

    # gpu LB sweep & boundaries
    generate_sweep(ctx, 'UniformGridGPU_LbKernel', update_rule,
                   field_swaps=[('pdfs', 'pdfs_tmp')],
                   inner_outer_split=True, target='gpu', gpu_indexing_params=sweep_params,
                   varying_parameters=vp)

    generate_boundary(ctx, 'UniformGridGPU_NoSlip', NoSlip(), lb_method, target='gpu')
    generate_boundary(ctx, 'UniformGridGPU_UBB', UBB([0.05, 0, 0]), lb_method, target='gpu')

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=1.0)
    getter_assignments = macroscopic_values_getter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=None)
    generate_sweep(ctx, 'UniformGridGPU_MacroSetter', setter_assignments)
    generate_sweep(ctx, 'UniformGridGPU_MacroGetter', getter_assignments)

    # communication
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_PackInfo', update_rule, target='gpu')

    infoHeaderParams = {
        'stencil': stencil_str,
        'q': q,
        'configName': ctx.config,
        'cse_global': int(options['optimization']['cse_global']),
        'cse_pdfs': int(options['optimization']['cse_pdfs']),
    }
    ctx.write_file("UniformGridGPU_Defines.h", info_header.format(**infoHeaderParams))
