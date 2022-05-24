import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.boundaries import NoSlip, UBB
from lbmpy import LBMConfig, LBMOptimisation, Stencil, Method, LBStencil
from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import RefinementScaling, generate_boundary, generate_lattice_model

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    omega, omega_free = sp.symbols("omega, omega_free")
    force_field, vel_field, omega_out = ps.fields(f"force(3), velocity(3), omega_out: {data_type}[3D]", layout='fzyx')

    stencil = LBStencil(Stencil.D3Q19)
    lbm_config = LBMConfig(stencil=stencil, method=Method.MRT, entropic=True, zero_centered=False,
                           compressible=True, omega_output_field=omega_out,
                           force=force_field.center_vector, output={'velocity': vel_field},
                           relaxation_rates=[omega, omega, omega_free, omega_free, omega_free, omega_free])

    lbm_opt = LBMOptimisation(cse_global=True)

    # the collision rule of the LB method where the some advanced features
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    # the refinement scaling object describes how certain parameters are scaled across grid scales
    # there are two default scaling behaviors available for relaxation rates and forces:
    scaling = RefinementScaling()
    scaling.add_standard_relaxation_rate_scaling(omega)
    scaling.add_force_scaling(force_field)

    # generate lattice model and (optionally) boundary conditions
    # for CPU simulations waLBerla's internal boundary handling can be used as well
    # If the field layout 'fzyx' is chosen vectorisation is usually possible. The default layout 'zyxf' allows only
    # for vectorisation on AVX512 due to scatter and gather intrinsics
    generate_lattice_model(ctx, 'LbCodeGenerationExample_LatticeModel', collision_rule,
                           field_layout='fzyx', refinement_scaling=scaling)
    generate_boundary(ctx, 'LbCodeGenerationExample_UBB', UBB([0.05, 0, 0], data_type=data_type), collision_rule.method)
    generate_boundary(ctx, 'LbCodeGenerationExample_NoSlip', NoSlip(), collision_rule.method)
    generate_info_header(ctx, 'LbCodeGenerationExample')
