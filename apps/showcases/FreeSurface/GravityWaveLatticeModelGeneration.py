import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, create_lb_collision_rule
from lbmpy.enums import ForceModel, Method, Stencil
from lbmpy.stencils import LBStencil

from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model


with CodeGeneration() as ctx:
    # general parameters
    layout = 'fzyx'
    data_type = "float64" if ctx.double_accuracy else "float32"

    stencil = LBStencil(Stencil.D3Q19)
    omega = sp.Symbol('omega')
    force_field = ps.fields(f"force(3): {data_type}[3D]", layout='fzyx')

    # method definition
    lbm_config = LBMConfig(stencil=stencil,
                           method=Method.SRT,
                           relaxation_rate=omega,
                           compressible=True,
                           force=force_field,
                           force_model=ForceModel.GUO,
                           zero_centered=False,
                           streaming_pattern='pull')  # free surface implementation only works with pull pattern

    # optimizations to be used by the code generator
    lbm_opt = LBMOptimisation(cse_global=True,
                              field_layout=layout)

    collision_rule = create_lb_collision_rule(lbm_config=lbm_config,
                                              lbm_optimisation=lbm_opt)

    generate_lattice_model(ctx, "GravityWaveLatticeModel", collision_rule, field_layout=layout)
