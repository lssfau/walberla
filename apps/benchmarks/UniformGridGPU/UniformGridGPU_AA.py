import sympy as sp
import numpy as np
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.fieldaccess import AAEvenTimeStepAccessor, AAOddTimeStepAccessor
from pystencils_walberla import generate_pack_info_from_kernel
from pystencils_walberla import CodeGeneration, generate_sweep
from pystencils.data_types import TypedSymbol
from pystencils.fast_approximation import insert_fast_sqrts, insert_fast_divisions
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter, macroscopic_values_setter

omega = sp.symbols("omega")
omega_free = sp.Symbol("omega_free")
compile_time_block_size = False

if compile_time_block_size:
    sweep_block_size = (128, 1, 1)
else:
    sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                        TypedSymbol("cudaBlockSize1", np.int32),
                        1)

sweep_params = {'block_size': sweep_block_size}

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
        'relaxation_rate': omega,
    },
    'mrt': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'relaxation_rates': [omega, 1.3, 1.4, omega, 1.2, 1.1],
    },
    'entropic': {
        'method': 'mrt',
        'stencil': 'D3Q19',
        'compressible': True,
        'relaxation_rates': [omega, omega, omega_free, omega_free, omega_free],
        'entropic': True,
    },
    'smagorinsky': {
        'method': 'srt',
        'stencil': 'D3Q19',
        'smagorinsky': True,
        'relaxation_rate': omega,
    }
}


info_header = """
#include "stencil/D3Q{q}.h"\nusing Stencil_T = walberla::stencil::D3Q{q};
const char * infoStencil = "{stencil}";
const char * infoConfigName = "{configName}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""


with CodeGeneration() as ctx:
    accessors = {
        'Even': AAEvenTimeStepAccessor(),
        'Odd': AAOddTimeStepAccessor()
    }

    common_options = {
        'field_name': 'pdfs',
        'optimization': {'cse_global': True,
                         'cse_pdfs': False,
                         'field_layout': 'fzyx',
                         }
    }
    options = options_dict.get(ctx.config, options_dict['srt'])
    options.update(common_options)

    stencil_str = options['stencil']
    q = int(stencil_str[stencil_str.find('Q') + 1:])
    pdfs, velocity_field = ps.fields("pdfs({q}), velocity(3) : double[3D]".format(q=q), layout='fzyx')
    options['optimization']['symbolic_field'] = pdfs

    vp = [
        ('int32_t', 'cudaBlockSize0'),
        ('int32_t', 'cudaBlockSize1')
    ]
    lb_method = create_lb_method(**options)

    # Kernels
    options_without_opt = options.copy()
    del options_without_opt['optimization']
    update_rules = {}
    for name, accessor in accessors.items():
        update_rule = create_lb_update_rule(lb_method=lb_method, kernel_type=accessor, **options)
        update_rule = insert_fast_divisions(update_rule)
        update_rule = insert_fast_sqrts(update_rule)
        update_rules[name] = update_rule
        generate_sweep(ctx, 'UniformGridGPU_AA_LbKernel' + name, update_rule,
                       inner_outer_split=True, target='gpu', gpu_indexing_params=sweep_params,
                       varying_parameters=vp)

    # getter & setter
    setter_assignments = macroscopic_values_setter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=1)
    getter_assignments = macroscopic_values_getter(lb_method, velocity=velocity_field.center_vector,
                                                   pdfs=pdfs.center_vector, density=None)
    generate_sweep(ctx, 'UniformGridGPU_AA_MacroSetter', setter_assignments)
    generate_sweep(ctx, 'UniformGridGPU_AA_MacroGetter', getter_assignments)

    # communication
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_AA_PackInfoPull', update_rules['Odd'],
                                   kind='pull', target='gpu')
    generate_pack_info_from_kernel(ctx, 'UniformGridGPU_AA_PackInfoPush', update_rules['Odd'],
                                   kind='push', target='gpu')

    infoHeaderParams = {
        'stencil': stencil_str,
        'q': q,
        'configName': ctx.config,
        'cse_global': int(options['optimization']['cse_global']),
        'cse_pdfs': int(options['optimization']['cse_pdfs']),
    }
    ctx.write_file("UniformGridGPU_AA_Defines.h", info_header.format(**infoHeaderParams))
