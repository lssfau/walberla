import sympy as sp
import numpy as np
import pystencils as ps

from pystencils.data_types import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions

from lbmpy.advanced_streaming import Timestep, is_inplace
from lbmpy.advanced_streaming.utility import streaming_patterns
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.stencils import get_stencil

from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep
from lbmpy_walberla import generate_alternating_lbm_sweep, generate_lb_pack_info, generate_alternating_lbm_boundary

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
compile_time_block_size = False

if compile_time_block_size:
    sweep_block_size = (128, 1, 1)
else:
    sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                        TypedSymbol("cudaBlockSize1", np.int32),
                        TypedSymbol("cudaBlockSize2", np.int32))

gpu_indexing_params = {'block_size': sweep_block_size}

options_dict = {
    'srt': {
        'method': 'srt',
        'relaxation_rate': omega,
        'compressible': False,
    },
    'trt': {
        'method': 'trt',
        'relaxation_rate': omega,
    },
    'mrt': {
        'method': 'mrt',
        'relaxation_rates': [omega, 1, 1, 1, 1, 1, 1],
    },
    'mrt-overrelax': {
        'method': 'mrt',
        'relaxation_rates': [omega, 1.3, 1.4, omega, 1.2, 1.1],
    },
    'cumulant': {
        'method': 'cumulant',
        'relaxation_rate': omega,
        'compressible': True,
    },
    'cumulant-overrelax': {
        'method': 'cumulant',
        'relaxation_rates': [omega] + [1 + x * 1e-2 for x in range(1, 11)],
        'compressible': True,
    },
    'entropic': {
        'method': 'mrt',
        'compressible': True,
        'relaxation_rates': [omega, omega, omega_free, omega_free, omega_free],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'smagorinsky': True,
        'relaxation_rate': omega,
    }
}

info_header = """
const char * infoStencil = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

# DEFAULTS
optimize = True

with CodeGeneration() as ctx:
    config_tokens = ctx.config.split('_')

    assert len(config_tokens) >= 3
    stencil_str = config_tokens[0]
    streaming_pattern = config_tokens[1]
    collision_setup = config_tokens[2]

    if len(config_tokens) >= 4:
        optimize = (config_tokens[3] != 'noopt')

    stencil = get_stencil(stencil_str)
    assert streaming_pattern in streaming_patterns, f"Invalid streaming pattern: {streaming_pattern}"

    options = options_dict[collision_setup]

    q = len(stencil)
    dim = len(stencil[0])
    assert dim == 3, "This app supports only three-dimensional stencils"
    pdfs, pdfs_tmp, velocity_field = ps.fields(f"pdfs({q}), pdfs_tmp({q}), velocity(3) : double[3D]", layout='fzyx')

    common_options = {
        'stencil': stencil,
        'field_name': pdfs.name,
        'optimization': {
            'target': 'gpu',
            'cse_global': True,
            'cse_pdfs': False,
            'symbolic_field': pdfs,
            'field_layout': 'fzyx',
            'gpu_indexing_params': gpu_indexing_params,
        }
    }

    options.update(common_options)

    if not is_inplace(streaming_pattern):
        options['optimization']['symbolic_temporary_field'] = pdfs_tmp
        field_swaps = [(pdfs, pdfs_tmp)]
    else:
        field_swaps = []

    vp = [
        ('int32_t', 'cudaBlockSize0'),
        ('int32_t', 'cudaBlockSize1'),
        ('int32_t', 'cudaBlockSize2')
    ]

    # LB Sweep
    collision_rule = create_lb_collision_rule(**options)

    if optimize:
        collision_rule = insert_fast_divisions(collision_rule)
        collision_rule = insert_fast_sqrts(collision_rule)

    lb_method = collision_rule.method

    generate_alternating_lbm_sweep(ctx, 'UniformGridGPU_LbKernel', collision_rule, streaming_pattern,
                                   optimization=options['optimization'],
                                   inner_outer_split=True, varying_parameters=vp, field_swaps=field_swaps)

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, density=1.0, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs,
                                                   streaming_pattern=streaming_pattern,
                                                   previous_timestep=Timestep.EVEN)
    generate_sweep(ctx, 'UniformGridGPU_MacroSetter', setter_assignments, target='gpu')

    # Boundaries
    noslip = NoSlip()
    ubb = UBB((0.05, 0, 0))

    generate_alternating_lbm_boundary(ctx, 'UniformGridGPU_NoSlip', noslip, lb_method, field_name=pdfs.name,
                                      streaming_pattern=streaming_pattern, target='gpu')
    generate_alternating_lbm_boundary(ctx, 'UniformGridGPU_UBB', ubb, lb_method, field_name=pdfs.name,
                                      streaming_pattern=streaming_pattern, target='gpu')

    # communication
    generate_lb_pack_info(ctx, 'UniformGridGPU_PackInfo', stencil, pdfs,
                          streaming_pattern=streaming_pattern, target='gpu',
                          always_generate_separate_classes=True)

    infoHeaderParams = {
        'stencil': stencil_str,
        'streaming_pattern': streaming_pattern,
        'collision_setup': collision_setup,
        'cse_global': int(options['optimization']['cse_global']),
        'cse_pdfs': int(options['optimization']['cse_pdfs']),
    }

    stencil_typedefs = {'Stencil_T': stencil,
                        'CommunicationStencil_T': stencil}
    field_typedefs = {'PdfField_T': pdfs,
                      'VelocityField_T': velocity_field}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'UniformGridGPU_InfoHeader',
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
