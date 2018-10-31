import sympy as sp
from lbmpy_walberla import generate_lattice_model_files
from lbmpy.creationfunctions import create_lb_update_rule
from pystencils_walberla.sweep import Sweep
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import create_lb_method
from lbmpy_walberla.boundary import create_boundary_class
from pystencils_walberla.cmake_integration import codegen


dtype = 'float64'

# LB options
options = {
    'method': 'srt',
    'stencil': 'D3Q19',
    'relaxation_rate': sp.Symbol("omega"),
    'field_name': 'pdfs',
    'compressible': False,
    'temporary_field_name': 'pdfs_tmp',
    'optimization': {'cse_global': True,
                     'cse_pdfs': True,
                     'double_precision': dtype == 'float64'}
}

# GPU optimization options
inner_opt = {'gpu_indexing_params': {'block_size': (128, 1, 1)},  'data_type': dtype}
outer_opt = {'gpu_indexing_params': {'block_size': (32, 32, 32)}, 'data_type': dtype}


def lb_assignments():
    ur = create_lb_update_rule(**options)
    return ur.all_assignments


def genBoundary():
    boundary = UBB([0.05, 0, 0], dim=3, name="UniformGridGPU_UBB")
    return create_boundary_class(boundary, create_lb_method(**options), target='gpu')


def genNoSlip():
    boundary = NoSlip(name='UniformGridGPU_NoSlip')
    return create_boundary_class(boundary, create_lb_method(**options), target='gpu')


generate_lattice_model_files(class_name='UniformGridGPU_LatticeModel', **options)

Sweep.generate_inner_outer_kernel('UniformGridGPU_LbKernel',
                                  lambda: create_lb_update_rule(**options).all_assignments,
                                  target='gpu',
                                  temporary_fields=['pdfs_tmp'],
                                  field_swaps=[('pdfs', 'pdfs_tmp')],
                                  optimization=inner_opt,
                                  outer_optimization=outer_opt)

Sweep.generate_pack_info('UniformGridGPU_PackInfo', lb_assignments, target='gpu')

codegen.register(['UniformGridGPU_UBB.h', 'UniformGridGPU_UBB.cu'], genBoundary)
codegen.register(['UniformGridGPU_NoSlip.h', 'UniformGridGPU_NoSlip.cu'], genNoSlip)
