import sympy as sp
import numpy as np

import pystencils as ps
from pystencils.typing import TypedSymbol

from lbmpy.advanced_streaming.utility import get_timesteps
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_method, create_lb_collision_rule
from lbmpy import LBMConfig, LBMOptimisation, Stencil, Method, LBStencil, SubgridScaleModel

from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
compile_time_block_size = False
max_threads = 256

sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                    TypedSymbol("cudaBlockSize1", np.int32),
                    TypedSymbol("cudaBlockSize2", np.int32))

gpu_indexing_params = {'block_size': sweep_block_size}

info_header = """
const char * infoStencil = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"

    streaming_pattern = 'esopull'
    timesteps = get_timesteps(streaming_pattern)
    stencil = LBStencil(Stencil.D3Q19)
    method_enum = Method.CUMULANT

    fourth_order_correction = 0.01 if method_enum == Method.CUMULANT and stencil.Q == 27 else False
    collision_setup = "cumulant-K17" if fourth_order_correction else method_enum.name.lower()

    assert stencil.D == 3, "This application supports only three-dimensional stencils"
    pdfs, pdfs_tmp = ps.fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {field_type}[3D]", layout='fzyx')
    density_field, velocity_field = ps.fields(f"density, velocity(3) : {field_type}[3D]", layout='fzyx')
    macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

    lbm_config = LBMConfig(stencil=stencil, method=method_enum, relaxation_rate=omega, compressible=True,
                           fourth_order_correction=fourth_order_correction,
                           streaming_pattern=streaming_pattern)
    lbm_opt = LBMOptimisation(cse_global=False, field_layout='fzyx')

    method = create_lb_method(lbm_config=lbm_config)
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip())
    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB([0.05, 0, 0], data_type=field_type))

    generate_lbm_package(ctx, name="NonUniformGridGPU",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=True, boundaries=[no_slip, ubb],
                         macroscopic_fields=macroscopic_fields,
                         target=ps.Target.GPU, gpu_indexing_params=gpu_indexing_params,
                         max_threads=max_threads)

    infoHeaderParams = {
        'stencil': stencil.name.lower(),
        'streaming_pattern': streaming_pattern,
        'collision_setup': collision_setup,
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }

    field_typedefs = {'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    generate_info_header(ctx, 'NonUniformGridGPUInfoHeader',
                         field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
