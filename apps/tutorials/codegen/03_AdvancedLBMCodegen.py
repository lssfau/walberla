import sympy as sp
import pystencils as ps

from lbmpy.creationfunctions import create_lb_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel
from lbmpy_walberla import generate_boundary


#   ========================
#      General Parameters
#   ========================

STENCIL = 'D2Q9'
OMEGA = sp.Symbol('omega')
LAYOUT = 'fzyx'

#   PDF Fields
pdfs, pdfs_tmp = ps.fields('pdfs(9), pdfs_tmp(9): [2D]', layout=LAYOUT)

#   Velocity Output Field
velocity = ps.fields("velocity(2): [2D]", layout=LAYOUT)
OUTPUT = {'velocity': velocity}

#   Optimization
OPT = {'target': 'cpu',
       'cse_global': True,
       'symbolic_field': pdfs,
       'symbolic_temporary_field': pdfs_tmp,
       'field_layout': LAYOUT}


#   ==================
#      Method Setup
#   ==================

lbm__params = {'stencil': STENCIL,
               'method': 'mrt_raw',
               'relaxation_rates': [0, 0, 0, OMEGA, OMEGA, OMEGA, 1, 1, 1],
               'cumulant': True,
               'compressible': True}

lbm_update_rule = create_lb_update_rule(optimization=OPT,
                                        output=OUTPUT,
                                        **lbm__params)

lbm_method = lbm_update_rule.method

#   ========================
#      PDF Initialization
#   ========================

initial_rho = sp.Symbol('rho_0')

pdfs_setter = macroscopic_values_setter(lbm_method, 
                                        initial_rho, 
                                        velocity.center_vector, 
                                        pdfs.center_vector)

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:

    #   LBM Sweep
    generate_sweep(ctx, "CumulantMRTSweep", lbm_update_rule, field_swaps=[(pdfs, pdfs_tmp)])

    #   Pack Info
    generate_pack_info_from_kernel(ctx, "CumulantMRTPackInfo", lbm_update_rule)

    #   Macroscopic Values Setter
    generate_sweep(ctx, "DensityAndVelocityFieldSetter", pdfs_setter)

    #   NoSlip Boundary
    generate_boundary(ctx, "CumulantMRTNoSlip", NoSlip(), lbm_method)
