import sympy as sp
import numpy as np
import pystencils as ps

from dataclasses import replace

from pystencils.typing import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.advanced_streaming import Timestep, is_inplace
from lbmpy.advanced_streaming.utility import streaming_patterns
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.moments import get_default_moment_set_for_stencil
from lbmpy.updatekernels import create_stream_only_kernel
from lbmpy.fieldaccess import *

from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep
from lbmpy_walberla import generate_alternating_lbm_sweep, generate_lb_pack_info, generate_alternating_lbm_boundary

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
compile_time_block_size = False
max_threads = None

if compile_time_block_size:
    sweep_block_size = (128, 1, 1)
else:
    sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                        TypedSymbol("cudaBlockSize1", np.int32),
                        TypedSymbol("cudaBlockSize2", np.int32))

gpu_indexing_params = {'block_size': sweep_block_size}

options_dict = {
    'srt': {
        'method': Method.SRT,
        'relaxation_rate': omega,
        'compressible': False,
    },
    'trt': {
        'method': Method.TRT,
        'relaxation_rate': omega,
        'compressible': False,
    },
    'mrt': {
        'method': Method.MRT,
        'relaxation_rates': [omega, 1, 1, 1, 1, 1, 1],
        'compressible': False,
    },
    'mrt-overrelax': {
        'method': Method.MRT,
        'relaxation_rates': [omega] + [1 + x * 1e-2 for x in range(1, 11)],
        'compressible': False,
    },
    'central': {
        'method': Method.CENTRAL_MOMENT,
        'relaxation_rate': omega,
        'compressible': True,
    },
    'central-overrelax': {
        'method': Method.CENTRAL_MOMENT,
        'relaxation_rates': [omega] + [1 + x * 1e-2 for x in range(1, 11)],
        'compressible': True,
    },
    'cumulant': {
        'method': Method.MONOMIAL_CUMULANT,
        'relaxation_rate': omega,
        'compressible': True,
    },
    'cumulant-overrelax': {
        'method': Method.MONOMIAL_CUMULANT,
        'relaxation_rates': [omega] + [1 + x * 1e-2 for x in range(1, 18)],
        'compressible': True,
    },
    'entropic': {
        'method': Method.TRT_KBC_N4,
        'compressible': True,
        'zero_centered': False,
        'relaxation_rates': [omega, omega_free],
        'entropic': True,
        'entropic_newton_iterations': False
    },
    'smagorinsky': {
        'method': Method.SRT,
        'smagorinsky': False,
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
    field_type = "float64" if ctx.double_accuracy else "float32"
    config_tokens = ctx.config.split('_')

    assert len(config_tokens) >= 3
    stencil_str = config_tokens[0]
    streaming_pattern = config_tokens[1]
    collision_setup = config_tokens[2]

    if len(config_tokens) >= 4:
        optimize = (config_tokens[3] != 'noopt')

    if stencil_str == "d3q27":
        stencil = LBStencil(Stencil.D3Q27)
    elif stencil_str == "d3q19":
        stencil = LBStencil(Stencil.D3Q19)
    else:
        raise ValueError("Only D3Q27 and D3Q19 stencil are supported at the moment")

    assert streaming_pattern in streaming_patterns, f"Invalid streaming pattern: {streaming_pattern}"

    options = options_dict[collision_setup]

    q = stencil.Q
    dim = stencil.D
    assert dim == 3, "This app supports only three-dimensional stencils"
    pdfs, pdfs_tmp, velocity_field = ps.fields(f"pdfs({q}), pdfs_tmp({q}), velocity(3) : {field_type}[3D]",
                                               layout='fzyx')

    lbm_config = LBMConfig(stencil=stencil, field_name=pdfs.name, streaming_pattern=streaming_pattern, **options)
    lbm_opt = LBMOptimisation(cse_global=True, cse_pdfs=False, symbolic_field=pdfs, field_layout='fzyx')

    if lbm_config.method == Method.CENTRAL_MOMENT:
        lbm_config = replace(lbm_config, nested_moments=get_default_moment_set_for_stencil(stencil))

    if not is_inplace(streaming_pattern):
        lbm_opt = replace(lbm_opt, symbolic_temporary_field=pdfs_tmp)
        field_swaps = [(pdfs, pdfs_tmp)]
    else:
        field_swaps = []

    vp = [
        ('int32_t', 'cudaBlockSize0'),
        ('int32_t', 'cudaBlockSize1'),
        ('int32_t', 'cudaBlockSize2')
    ]

    # Sweep for Stream only. This is for benchmarking an empty streaming pattern without LBM.
    # is_inplace is set to False to ensure that the streaming is done with src and dst field.
    # If this is not the case the compiler might simplify the streaming in a way that benchmarking makes no sense.
    accessor = CollideOnlyInplaceAccessor()
    accessor.is_inplace = False
    field_swaps_stream_only = [(pdfs, pdfs_tmp)]
    stream_only_kernel = create_stream_only_kernel(stencil, pdfs, pdfs_tmp, accessor=accessor)

    # LB Sweep
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    if optimize:
        collision_rule = insert_fast_divisions(collision_rule)
        collision_rule = insert_fast_sqrts(collision_rule)

    lb_method = collision_rule.method

    generate_alternating_lbm_sweep(ctx, 'UniformGridGPU_LbKernel', collision_rule, lbm_config=lbm_config,
                                   lbm_optimisation=lbm_opt, target=ps.Target.GPU,
                                   gpu_indexing_params=gpu_indexing_params,
                                   inner_outer_split=True, varying_parameters=vp, field_swaps=field_swaps,
                                   max_threads=max_threads)

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, density=1.0, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs,
                                                   streaming_pattern=streaming_pattern,
                                                   previous_timestep=Timestep.EVEN)
    generate_sweep(ctx, 'UniformGridGPU_MacroSetter', setter_assignments, target=ps.Target.GPU, max_threads=max_threads)

    # Stream only kernel
    generate_sweep(ctx, 'UniformGridGPU_StreamOnlyKernel', stream_only_kernel, field_swaps=field_swaps_stream_only,
                   gpu_indexing_params=gpu_indexing_params, varying_parameters=vp, target=ps.Target.GPU,
                   max_threads=max_threads)

    # Boundaries
    noslip = NoSlip()
    ubb = UBB((0.05, 0, 0), data_type=field_type)

    generate_alternating_lbm_boundary(ctx, 'UniformGridGPU_NoSlip', noslip, lb_method, field_name=pdfs.name,
                                      streaming_pattern=streaming_pattern, target=ps.Target.GPU)
    generate_alternating_lbm_boundary(ctx, 'UniformGridGPU_UBB', ubb, lb_method, field_name=pdfs.name,
                                      streaming_pattern=streaming_pattern, target=ps.Target.GPU)

    # communication
    generate_lb_pack_info(ctx, 'UniformGridGPU_PackInfo', stencil, pdfs,
                          streaming_pattern=streaming_pattern, target=ps.Target.GPU,
                          always_generate_separate_classes=True)

    infoHeaderParams = {
        'stencil': stencil_str,
        'streaming_pattern': streaming_pattern,
        'collision_setup': collision_setup,
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }

    stencil_typedefs = {'Stencil_T': stencil,
                        'CommunicationStencil_T': stencil}
    field_typedefs = {'PdfField_T': pdfs,
                      'VelocityField_T': velocity_field}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'UniformGridGPU_InfoHeader',
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
