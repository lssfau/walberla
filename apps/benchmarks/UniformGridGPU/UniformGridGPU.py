import sympy as sp
import numpy as np
import pystencils as ps

from dataclasses import replace

from pystencils import Assignment
from pystencils.typing import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, SubgridScaleModel
from lbmpy.advanced_streaming import is_inplace
from lbmpy.advanced_streaming.utility import streaming_patterns
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.moments import get_default_moment_set_for_stencil

from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
compile_time_block_size = False
max_threads = 256

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
    'cumulant-K17': {
        'method': Method.CUMULANT,
        'relaxation_rate': omega,
        'compressible': True,
        'fourth_order_correction': 0.01
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
        'subgrid_scale_model': SubgridScaleModel.SMAGORINSKY,
        'relaxation_rate': omega,
    },
    'qr': {
        'method': Method.SRT,
        'subgrid_scale_model': SubgridScaleModel.QR,
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
    pdf_data_type = "float64"
    field_data_type = "float64"
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

    assert stencil.D == 3, "This application supports only three-dimensional stencils"
    pdfs, pdfs_tmp = ps.fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {pdf_data_type}[3D]", layout='fzyx')
    density_field, velocity_field = ps.fields(f"density, velocity(3) : {field_data_type}[3D]", layout='fzyx')
    macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

    lbm_config = LBMConfig(stencil=stencil, field_name=pdfs.name, streaming_pattern=streaming_pattern, **options)
    lbm_opt = LBMOptimisation(cse_global=True, cse_pdfs=False, symbolic_field=pdfs, field_layout='fzyx')

    if lbm_config.method == Method.CENTRAL_MOMENT:
        lbm_config = replace(lbm_config, nested_moments=get_default_moment_set_for_stencil(stencil))

    if not is_inplace(streaming_pattern):
        lbm_opt = replace(lbm_opt, symbolic_temporary_field=pdfs_tmp)
        field_swaps = [(pdfs, pdfs_tmp)]
    else:
        field_swaps = []

    # This is a microbenchmark for testing how fast Q PDFs can be updated per cell. To avoid optimisations from
    # the compiler the PDFs are shuffled inside a cell. Otherwise, for common streaming patterns compilers would
    # typically remove the copy of the center PDF which results in an overestimation of the maximum performance
    stream_only_kernel = []
    for i in range(stencil.Q):
        stream_only_kernel.append(Assignment(pdfs(i), pdfs((i + 3) % stencil.Q)))

    # LB Sweep
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    if optimize:
        collision_rule = insert_fast_divisions(collision_rule)
        collision_rule = insert_fast_sqrts(collision_rule)

    lb_method = collision_rule.method

    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip(), field_data_type=pdf_data_type)
    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB([0.05, 0, 0], data_type=field_data_type),
                                 field_data_type=pdf_data_type)

    generate_lbm_package(ctx, name="UniformGridGPU",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=False, boundaries=[no_slip, ubb],
                         macroscopic_fields=macroscopic_fields,
                         target=ps.Target.GPU, gpu_indexing_params=gpu_indexing_params,
                         data_type=field_data_type, pdfs_data_type=pdf_data_type,
                         max_threads=max_threads)

    # Stream only kernel
    generate_sweep(ctx, 'UniformGridGPU_StreamOnlyKernel', stream_only_kernel,
                   gpu_indexing_params={'block_size': (128, 1, 1)}, target=ps.Target.GPU,
                   max_threads=max_threads)

    infoHeaderParams = {
        'stencil': stencil_str,
        'streaming_pattern': streaming_pattern,
        'collision_setup': collision_setup,
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }

    field_typedefs = {'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'UniformGridGPU_InfoHeader',
                         field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
