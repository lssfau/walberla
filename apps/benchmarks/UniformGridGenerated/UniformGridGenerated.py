import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_update_rule, create_lb_collision_rule
from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel, generate_sweep,\
    generate_mpidtype_info_from_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter, macroscopic_values_setter
from lbmpy.fieldaccess import AAEvenTimeStepAccessor, AAOddTimeStepAccessor

omega = sp.symbols('omega')
omega_free = sp.Symbol('omega_free')
omega_fill = sp.symbols('omega_:10')

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
        'relaxation_rates': [0.0, omega, 1.3, 1.4, omega, 1.2, 1.1, 1.15, 1.234, 1.4235, 1.242, 1.2567, 0.9, 0.7],
    },
    'mrt_full': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [omega_fill[0], omega, omega_fill[1], omega_fill[2],
                             omega_fill[3], omega_fill[4], omega_fill[5]],
    },
    'entropic': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rates': [omega, omega, omega_free, omega_free, omega_free],
        'entropic': True,
    },
    'entropic_kbc_n4': {
        'method': 'trt-kbc-n4',
        'stencil': 'D3Q27',
        'compressible': True,
        'relaxation_rates': [omega, omega_free],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'smagorinsky': True,
        'relaxation_rate': omega,
    },
    'cumulant': {
        'method': 'cumulant',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rate': omega,
    },
}

info_header = """
#include "stencil/D3Q{q}.h"\nusing Stencil_T = walberla::stencil::D3Q{q};
const char * infoStencil = "{stencil}";
const char * infoConfigName = "{configName}";
const char * optimizationDict = "{optimizationDict}";
"""

with CodeGeneration() as ctx:
    common_options = {
        'field_name': 'pdfs',
        'temporary_field_name': 'pdfs_tmp',
    }
    opts = {
        'two_field_cse_pdfs': False,
        'two_field_cse_global': False,
        'two_field_split': True,
        'two_field_nt_stores': True,

        'aa_even_cse_pdfs': False,
        'aa_even_cse_global': False,
        'aa_even_split': False,
        'aa_even_nt_stores': False,

        'aa_odd_cse_pdfs': False,
        'aa_odd_cse_global': False,
        'aa_odd_split': True,
        'aa_odd_nt_stores': False,

        'compiled_in_boundaries': False,
    }
    config_name = ctx.config
    noopt = False
    d3q27 = False
    if config_name.endswith('_d3q27'):
        d3q27 = True
        config_name = config_name[:-len('_d3q27')]

    if config_name == '':
        config_name = 'trt'
    options = options_dict[config_name]
    options.update(common_options)
    options = options.copy()

    if d3q27:
        options['stencil'] = 'D3Q27'

    dtype_string = 'float64' if ctx.double_accuracy else 'float32'

    stencil_str = options['stencil']
    q = int(stencil_str[stencil_str.find('Q') + 1:])
    pdfs, velocity_field = ps.fields(f'pdfs({q}), velocity(3) : {dtype_string}[3D]', layout='fzyx')

    update_rule_two_field = create_lb_update_rule(optimization={'symbolic_field': pdfs,
                                                                'split': opts['two_field_split'],
                                                                'cse_global': opts['two_field_cse_global'],
                                                                'cse_pdfs': opts['two_field_cse_pdfs']}, **options)

    if opts['compiled_in_boundaries']:
        from lbmpy.boundaries import NoSlip, UBB
        from lbmpy.boundaries.boundaries_in_kernel import update_rule_with_push_boundaries
        from collections import OrderedDict
        boundaries = OrderedDict((
            ((1, 0, 0), NoSlip()),
            ((-1, 0, 0), NoSlip()),
            ((0, 1, 0), NoSlip()),
            ((0, -1, 0), NoSlip()),
            ((0, 0, 1), UBB([0.05, 0, 0])),
            ((0, 0, -1), NoSlip()),
        ))
        cr_even = create_lb_collision_rule(stencil='D3Q19', compressible=False,
                                           optimization={'cse_global': opts['aa_even_cse_global'],
                                                         'cse_pdfs': opts['aa_even_cse_pdfs']})
        cr_odd = create_lb_collision_rule(stencil='D3Q19', compressible=False,
                                          optimization={'cse_global': opts['aa_odd_cse_global'],
                                                        'cse_pdfs': opts['aa_odd_cse_pdfs']})
        update_rule_aa_even = update_rule_with_push_boundaries(cr_even, pdfs, boundaries,
                                                               AAEvenTimeStepAccessor, AAOddTimeStepAccessor.read)
        update_rule_aa_odd = update_rule_with_push_boundaries(cr_odd, pdfs, boundaries,
                                                              AAOddTimeStepAccessor, AAEvenTimeStepAccessor.read)
    else:
        update_rule_aa_even = create_lb_update_rule(kernel_type=AAEvenTimeStepAccessor(),
                                                    optimization={'symbolic_field': pdfs,
                                                                  'split': opts['aa_even_split'],
                                                                  'cse_global': opts['aa_even_cse_global'],
                                                                  'cse_pdfs': opts['aa_even_cse_pdfs']}, **options)
        update_rule_aa_odd = create_lb_update_rule(kernel_type=AAOddTimeStepAccessor(),
                                                   optimization={'symbolic_field': pdfs,
                                                                 'split': opts['aa_odd_split'],
                                                                 'cse_global': opts['aa_odd_cse_global'],
                                                                 'cse_pdfs': opts['aa_odd_cse_pdfs']}, **options)

    vec = {'assume_aligned': True, 'assume_inner_stride_one': True}

    # check if openmp is enabled in cmake
    if ctx.openmp:
        openmp_enabled = True
    else:
        openmp_enabled = False

    # Sweeps
    vec['nontemporal'] = opts['two_field_nt_stores']
    generate_sweep(ctx, 'GenLbKernel', update_rule_two_field, field_swaps=[('pdfs', 'pdfs_tmp')],
                   cpu_vectorize_info=vec)
    vec['nontemporal'] = opts['aa_even_nt_stores']
    generate_sweep(ctx, 'GenLbKernelAAEven', update_rule_aa_even, cpu_vectorize_info=vec,
                   cpu_openmp=openmp_enabled, ghost_layers=1)
    vec['nontemporal'] = opts['aa_odd_nt_stores']
    generate_sweep(ctx, 'GenLbKernelAAOdd', update_rule_aa_odd, cpu_vectorize_info=vec,
                   cpu_openmp=openmp_enabled, ghost_layers=1)

    setter_assignments = macroscopic_values_setter(update_rule_two_field.method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=1.0)
    getter_assignments = macroscopic_values_getter(update_rule_two_field.method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=None)
    generate_sweep(ctx, 'GenMacroSetter', setter_assignments, cpu_openmp=openmp_enabled)
    generate_sweep(ctx, 'GenMacroGetter', getter_assignments, cpu_openmp=openmp_enabled)

    # Communication
    generate_pack_info_from_kernel(ctx, 'GenPackInfo', update_rule_two_field,
                                   cpu_vectorize_info={'instruction_set': None}, cpu_openmp=False)
    generate_pack_info_from_kernel(ctx, 'GenPackInfoAAPull', update_rule_aa_odd, kind='pull',
                                   cpu_vectorize_info={'instruction_set': None}, cpu_openmp=False)
    generate_pack_info_from_kernel(ctx, 'GenPackInfoAAPush', update_rule_aa_odd, kind='push',
                                   cpu_vectorize_info={'instruction_set': None}, cpu_openmp=False)

    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfo', update_rule_two_field)
    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfoAAPull', update_rule_aa_odd, kind='pull')
    generate_mpidtype_info_from_kernel(ctx, 'GenMpiDtypeInfoAAPush', update_rule_aa_odd, kind='push')

    # Info Header
    infoHeaderParams = {
        'stencil': stencil_str,
        'q': q,
        'configName': ctx.config,
        'optimizationDict': str(opts),
    }
    ctx.write_file('GenDefines.h', info_header.format(**infoHeaderParams))
