import sympy as sp
from lbmpy_walberla import generate_lattice_model_files
from lbmpy.creationfunctions import create_lb_update_rule
from pystencils_walberla.sweep import Sweep

# LB options
options = {
    'method': 'srt',
    'stencil': 'D3Q19',
    'relaxation_rate': sp.Symbol("omega"),
    'field_name': 'pdfs',
    'compressible': False,
    'maxwellian_moments': False,
    'temporary_field_name': 'pdfs_tmp',
    'optimization': {'cse_global': False,
                     'cse_pdfs': False,
                     'double_precision': True}
}

# GPU optimization options
opt =       {'gpu_indexing_params': {'block_size': (128, 2, 1)},  'data_type': 'float64'}
outer_opt = {'gpu_indexing_params': {'block_size': (32, 32, 32)}, 'data_type': 'float64'}


def lb_assignments():
    ur = create_lb_update_rule(**options)
    return ur.all_assignments


generate_lattice_model_files(class_name='EquivalenceTest_LatticeModel', **options)

Sweep.generate_inner_outer_kernel('EquivalenceTest_GPUKernel',
                                  lambda: create_lb_update_rule(**options).all_assignments,
                                  target='gpu',
                                  temporary_fields=['pdfs_tmp'],
                                  field_swaps=[('pdfs', 'pdfs_tmp')],
                                  optimization=opt,
                                  outer_optimization=outer_opt)

Sweep.generate_pack_info('EquivalenceTest_GPUPackInfo', lb_assignments, target='gpu')
