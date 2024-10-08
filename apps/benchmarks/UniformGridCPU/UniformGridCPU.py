from dataclasses import replace

import sympy as sp
import pystencils as ps

from lbmpy.advanced_streaming import is_inplace
from lbmpy.advanced_streaming.utility import streaming_patterns
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, LBStencil, create_lb_collision_rule
from lbmpy.enums import Method, Stencil, SubgridScaleModel
from lbmpy.moments import get_default_moment_set_for_stencil

from pystencils_walberla import CodeGeneration, generate_info_header, generate_sweep
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

omega = sp.symbols('omega')
omega_free = sp.Symbol('omega_free')

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
const bool vectorised = {vec};
const bool nontemporal = {nt_stores};
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    openmp = True if ctx.openmp else False
    field_type = "float64" if ctx.double_accuracy else "float32"
    # This base pointer specification causes introduces temporary pointers in the outer loop such that the inner loop
    # only contains aligned memory addresses. Doing so NT Stores are much more effective which causes great perfomance
    # gains especially for the pull scheme on skylake architectures
    base_pointer_spec = None  # [['spatialInner0'], ['spatialInner1']]
    # cpu_vec = {"instruction_set": "best", "nontemporal": False,
    #            "assume_aligned": True, 'assume_sufficient_line_padding': True}

    cpu_vec = {"instruction_set": None}
    nt_stores = False

    config_tokens = ctx.config.split('_')

    assert len(config_tokens) >= 3
    stencil_str = config_tokens[0]
    streaming_pattern = config_tokens[1]
    collision_setup = config_tokens[2]

    if stencil_str == "d3q27":
        stencil = LBStencil(Stencil.D3Q27)
    elif stencil_str == "d3q19":
        stencil = LBStencil(Stencil.D3Q19)
    else:
        raise ValueError("Only D3Q27 and D3Q19 stencil are supported at the moment")

    assert streaming_pattern in streaming_patterns, f"Invalid streaming pattern: {streaming_pattern}"
    options = options_dict[collision_setup]

    assert stencil.D == 3, "This application supports only three-dimensional stencils"
    pdfs, pdfs_tmp = ps.fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {field_type}[3D]", layout='fzyx')
    density_field, velocity_field = ps.fields(f"density, velocity(3) : {field_type}[3D]", layout='fzyx')
    macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

    lbm_config = LBMConfig(stencil=stencil, field_name=pdfs.name, streaming_pattern=streaming_pattern, **options)
    lbm_opt = LBMOptimisation(cse_global=True, cse_pdfs=False, symbolic_field=pdfs, field_layout='fzyx')

    # This creates a simplified version of the central moment collision operator where the bulk and shear viscosity is
    # not seperated. This is done to get a fair comparison with the monomial cumulants.
    if lbm_config.method == Method.CENTRAL_MOMENT:
        lbm_config = replace(lbm_config, nested_moments=get_default_moment_set_for_stencil(stencil))

    if not is_inplace(streaming_pattern):
        lbm_opt = replace(lbm_opt, symbolic_temporary_field=pdfs_tmp)

    # This is a microbenchmark for testing how fast Q PDFs can be updated per cell. To avoid optimisations from
    # the compiler the PDFs are shuffled inside a cell. Otherwise, for common streaming patterns compilers would
    # typically remove the copy of the center PDF which results in an overestimation of the maximum performance
    stream_only_kernel = []
    for i in range(stencil.Q):
        stream_only_kernel.append(ps.Assignment(pdfs(i), pdfs((i + 3) % stencil.Q)))

    # LB Sweep
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip())
    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB([0.05, 0, 0], data_type=field_type))

    generate_lbm_package(ctx, name="UniformGridCPU",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=False, boundaries=[no_slip, ubb],
                         macroscopic_fields=macroscopic_fields,
                         cpu_openmp=openmp, cpu_vectorize_info=cpu_vec,
                         base_pointer_specification=base_pointer_spec)

    # Stream only kernel
    cpu_vec_stream = None
    if ctx.optimize_for_localhost:
        cpu_vec_stream = {"instruction_set": "best", "nontemporal": True,
                          "assume_aligned": True, 'assume_sufficient_line_padding': True,
                          "assume_inner_stride_one": True}

    generate_sweep(ctx, 'UniformGridCPU_StreamOnlyKernel', stream_only_kernel,
                   target=ps.Target.CPU, cpu_openmp=openmp,
                   cpu_vectorize_info=cpu_vec_stream, base_pointer_specification=[['spatialInner0'], ['spatialInner1']])

    infoHeaderParams = {
        'stencil': stencil_str,
        'streaming_pattern': streaming_pattern,
        'collision_setup': collision_setup,
        'vec': int(True if cpu_vec else False),
        'nt_stores': int(nt_stores),
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }

    field_typedefs = {'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'UniformGridCPU_InfoHeader',
                         field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
