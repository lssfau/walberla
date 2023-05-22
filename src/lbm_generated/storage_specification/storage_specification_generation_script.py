import sympy as sp

from pystencils import Target

from lbmpy.creationfunctions import create_lb_method
from lbmpy import LBMConfig, Stencil, Method, LBStencil
from pystencils_walberla import ManualCodeGenerationContext, generate_info_header
from lbmpy_walberla.storage_specification import generate_lbm_storage_specification


with ManualCodeGenerationContext(openmp=False, optimize_for_localhost=False,
                                 mpi=True, double_accuracy=True, cuda=False) as ctx:

    for stencil in [LBStencil(Stencil.D3Q19), LBStencil(Stencil.D3Q27)]:
        target = Target.GPU if ctx.cuda else Target.CPU
        data_type = "float64" if ctx.double_accuracy else "float32"

        method = Method.SRT
        relaxation_rate = sp.symbols("omega")
        streaming_pattern = 'pull'
        nonuniform = False

        lbm_config = LBMConfig(stencil=stencil, method=method, relaxation_rate=relaxation_rate,
                               streaming_pattern=streaming_pattern)

        lb_method = create_lb_method(lbm_config=lbm_config)

        storage_spec_name = f'{stencil.name}StorageSpecification'
        generate_lbm_storage_specification(ctx, storage_spec_name, lb_method, lbm_config,
                                           nonuniform=nonuniform, target=target, data_type=data_type)

        ctx.write_all_files()
