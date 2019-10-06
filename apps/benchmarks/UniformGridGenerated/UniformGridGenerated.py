import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_update_rule
from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel, generate_sweep, generate_mpidtype_info_from_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter, macroscopic_values_setter
from lbmpy.fieldaccess import AAEvenTimeStepAccessor, AAOddTimeStepAccessor

omega = sp.symbols("omega")
omega_fill = sp.symbols("omega_:10")

options_dict = {
    'srt': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'relaxation_rate': omega,
        'compressible': False,
    },
    'trt': {
        'method': 'trt',
        'stencil': 'D3Q19',
        'compressible': False,
        'relaxation_rate': omega,
    },
    'mrt': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [0, omega, 1.3, 1.4, omega, 1.2, 1.1, 1.15, 1.234, 1.4235, 1.242, 1.2567, 0.9, 0.7],
    },
    'mrt_full': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [omega_fill[0], omega, omega_fill[1], omega_fill[2], omega_fill[3], omega_fill[4], omega_fill[5]],
    },
    'mrt3': {
        'method': 'mrt3',
        'stencil': 'D3Q19',
        'relaxation_rates': [omega, 1.1, 1.2],
    },
    'entropic': {
        'method': 'mrt3',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rates': [omega, omega, sp.Symbol("omega_free")],
        'entropic': True,
    },
    'entropic_kbc_n4': {
        'method': 'trt-kbc-n4',
        'stencil': 'D3Q27',
        'compressible': True,
        'relaxation_rates': [omega, sp.Symbol("omega_free")],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'smagorinsky': True,
        'relaxation_rate': omega,
    },
    'cumulant': {
        'stencil': 'D3Q19',
        'compressible': True,
        'method': 'mrt',
        'cumulant': True,
        'relaxation_rates': [0, omega, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    },
}

info_header = """
#include "stencil/D3Q{q}.h"\nusing Stencil_T = walberla::stencil::D3Q{q}; 
const char * infoStencil = "{stencil}";
const char * infoConfigName = "{configName}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""


with CodeGeneration() as ctx:
    common_options = {
        'field_name': 'pdfs',
        'temporary_field_name': 'pdfs_tmp',
        'optimization': {'cse_global': False,
                         'cse_pdfs': False,
                         'split': False}
    }
    config_name = ctx.config
    noopt = False
    d3q27 = False
    if config_name.endswith("_noopt"):
        noopt = True
        config_name = config_name[:-len("_noopt")]
    if config_name.endswith("_d3q27"):
        d3q27 = True
        config_name = config_name[:-len("_d3q27")]

    if config_name == '':
        config_name = 'trt'
    options = options_dict[config_name]
    options.update(common_options)
    options = options.copy()

    if d3q27:
        options['stencil'] = 'D3Q27'

    stencil_str = options['stencil']
    q = int(stencil_str[stencil_str.find('Q')+1:])
    pdfs, velocity_field = ps.fields("pdfs({q}), velocity(3) : double[3D]".format(q=q), layout='fzyx')
    options['optimization']['symbolic_field'] = pdfs

    update_rule_two_field = create_lb_update_rule(**options)
    update_rule_aa_even = create_lb_update_rule(kernel_type=AAEvenTimeStepAccessor(), **options)
    options['optimization']['split'] = True
    update_rule_aa_odd = create_lb_update_rule(kernel_type=AAOddTimeStepAccessor(), **options)

    vec = {'nontemporal': False, 'assume_aligned': True, 'assume_inner_stride_one': True}

    # Sweeps
    generate_sweep(ctx, 'GenLbKernel', update_rule_two_field, field_swaps=[('pdfs', 'pdfs_tmp')])
    generate_sweep(ctx, 'GenLbKernelAAEven', update_rule_aa_even, cpu_vectorize_info={'assume_aligned': True},
                   cpu_openmp=True, ghost_layers=1)
    generate_sweep(ctx, 'GenLbKernelAAOdd', update_rule_aa_odd, cpu_vectorize_info={'assume_aligned': True},
                   cpu_openmp=True, ghost_layers=1)

    setter_assignments = macroscopic_values_setter(update_rule_two_field.method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=1)
    getter_assignments = macroscopic_values_getter(update_rule_two_field.method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=None)
    generate_sweep(ctx, 'GenMacroSetter', setter_assignments, cpu_openmp=True)
    generate_sweep(ctx, 'GenMacroGetter', getter_assignments, cpu_openmp=True)

    # Communication
    generate_pack_info_from_kernel(ctx, 'GenPackInfo', update_rule_two_field,
                                   cpu_vectorize_info={'instruction_set': None})
    generate_pack_info_from_kernel(ctx, 'GenPackInfoAAPull', update_rule_aa_odd, kind='pull',
                                   cpu_vectorize_info={'instruction_set': None})
    generate_pack_info_from_kernel(ctx, 'GenPackInfoAAPush', update_rule_aa_odd, kind='push',
                                   cpu_vectorize_info={'instruction_set': None})

    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfo', update_rule_two_field)
    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfoAAPull', update_rule_aa_odd, kind='pull')
    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfoAAPush', update_rule_aa_odd, kind='push')

    # Info Header
    infoHeaderParams = {
        'stencil': stencil_str,
        'q': q,
        'configName': ctx.config,
        'cse_global': int(options['optimization']['cse_global']),
        'cse_pdfs': int(options['optimization']['cse_pdfs']),
    }
    ctx.write_file("GenDefines.h", info_header.format(**infoHeaderParams))

