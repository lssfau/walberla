import sympy as sp
import pystencils as ps

from lbmpy import Stencil, LBStencil
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, Method, create_lb_collision_rule

from lbmpy_walberla import RefinementScaling, generate_lattice_model
from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel


with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q19)
    omega = sp.Symbol('omega')
    layout = 'fzyx'
    pdf_field = ps.fields(f"pdfs({stencil.Q}): {data_type}[3D]", layout=layout)

    lbm_config = LBMConfig(stencil=stencil, method=Method.CUMULANT, relaxation_rate=omega, compressible=True)
    lbm_opt = LBMOptimisation(cse_global=True, cse_pdfs=False, split=False,
                              symbolic_field=pdf_field, field_layout=layout, pre_simplification=True)

    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    refinement_scaling = RefinementScaling()
    refinement_scaling.add_standard_relaxation_rate_scaling(omega)

    generate_lattice_model(ctx, "LbCodeGen_LatticeModel", collision_rule, refinement_scaling=refinement_scaling,
                           field_layout=layout)

