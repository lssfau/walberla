import sympy as sp

from pystencils import Target
from pystencils import fields

from lbmpy.advanced_streaming.utility import get_timesteps
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_method, create_lb_collision_rule
from lbmpy import LBMConfig, LBMOptimisation, Stencil, Method, LBStencil
from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

import warnings

warnings.filterwarnings("ignore")
with CodeGeneration() as ctx:
    target = Target.CPU  # Target.GPU if ctx.cuda else Target.CPU
    data_type = "float64" if ctx.double_accuracy else "float32"

    streaming_pattern = 'esotwist'
    timesteps = get_timesteps(streaming_pattern)

    omega = sp.symbols("omega")

    stencil = LBStencil(Stencil.D3Q19)
    pdfs, vel_field = fields(f"pdfs({stencil.Q}), velocity({stencil.D}): {data_type}[{stencil.D}D]", layout='fzyx')

    macroscopic_fields = {'velocity': vel_field}

    lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=omega,
                           streaming_pattern=streaming_pattern)
    lbm_opt = LBMOptimisation(cse_global=False, field_layout='fzyx')

    method = create_lb_method(lbm_config=lbm_config)
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip())
    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB([0.05, 0, 0], data_type=data_type))

    generate_lbm_package(ctx, name="LBM",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=True, boundaries=[no_slip, ubb],
                         macroscopic_fields=macroscopic_fields)

    generate_info_header(ctx, 'Example_InfoHeader')
