import sympy as sp
import pystencils as ps

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil

from lbmpy.creationfunctions import create_lb_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel
from lbmpy_walberla import generate_boundary

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:
    #   ========================
    #      General Parameters
    #   ========================
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D2Q9)
    omega = sp.Symbol('omega')
    layout = 'fzyx'

    #   PDF Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[2D]', layout=layout)

    #   Velocity Output Field
    velocity = ps.fields(f"velocity({stencil.D}): {data_type}[2D]", layout=layout)
    output = {'velocity': velocity}

    # LBM Optimisation
    lbm_opt = LBMOptimisation(cse_global=True,
                              symbolic_field=pdfs,
                              symbolic_temporary_field=pdfs_tmp,
                              field_layout=layout)

    #   ==================
    #      Method Setup
    #   ==================

    lbm_config = LBMConfig(stencil=stencil,
                           method=Method.CUMULANT,
                           relaxation_rate=omega,
                           compressible=True,
                           output=output)

    lbm_update_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    lbm_method = lbm_update_rule.method

    #   ========================
    #      PDF Initialization
    #   ========================

    initial_rho = sp.Symbol('rho_0')

    pdfs_setter = macroscopic_values_setter(lbm_method,
                                            initial_rho,
                                            velocity.center_vector,
                                            pdfs.center_vector)

    target = ps.Target.GPU if ctx.gpu else ps.Target.CPU

    #   LBM Sweep
    generate_sweep(ctx, "CumulantMRTSweep", lbm_update_rule, field_swaps=[(pdfs, pdfs_tmp)], target=target)

    #   Pack Info
    generate_pack_info_from_kernel(ctx, "CumulantMRTPackInfo", lbm_update_rule, target=target)

    #   Macroscopic Values Setter
    generate_sweep(ctx, "InitialPDFsSetter", pdfs_setter, target=target)

    #   NoSlip Boundary
    generate_boundary(ctx, "CumulantMRTNoSlip", NoSlip(), lbm_method, target=target)
