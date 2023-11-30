import sympy as sp

import pystencils as ps

from lbmpy.advanced_streaming.utility import get_timesteps
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_method, create_lb_collision_rule
from lbmpy import LBMConfig, LBMOptimisation, Stencil, Method, LBStencil

from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")

info_header = """
const char * infoStencil = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"

    streaming_pattern = 'aa'
    timesteps = get_timesteps(streaming_pattern)
    stencil = LBStencil(Stencil.D3Q19)

    assert stencil.D == 3, "This application supports only three-dimensional stencils"
    pdfs, pdfs_tmp = ps.fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {field_type}[3D]", layout='fzyx')
    density_field, velocity_field = ps.fields(f"density, velocity(3) : {field_type}[3D]", layout='fzyx')
    macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

    lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=omega, compressible=True,
                           streaming_pattern=streaming_pattern)
    lbm_opt = LBMOptimisation(cse_global=False, field_layout="fzyx")

    method = create_lb_method(lbm_config=lbm_config)
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip())
    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB([0.05, 0, 0], data_type=field_type))

    generate_lbm_package(ctx, name="NonUniformGridCPU",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=True, boundaries=[no_slip, ubb],
                         macroscopic_fields=macroscopic_fields,
                         target=ps.Target.CPU)

    infoHeaderParams = {
        'stencil': stencil.name.lower(),
        'streaming_pattern': streaming_pattern,
        'collision_setup': lbm_config.method.name.lower(),
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }

    field_typedefs = {'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    generate_info_header(ctx, 'NonUniformGridCPUInfoHeader',
                         field_typedefs=field_typedefs,
                         additional_code=info_header.format(**infoHeaderParams))
