import sympy as sp

from pystencils import Target
from pystencils import fields

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy import LBMConfig, LBMOptimisation, Stencil, Method, LBStencil
from pystencils_walberla import ManualCodeGenerationContext, generate_info_header
from lbmpy_walberla import generate_lbm_sweep_collection


with ManualCodeGenerationContext(openmp=False, optimize_for_localhost=False,
                                 mpi=True, double_accuracy=True, cuda=False) as ctx:

    for stencil in [LBStencil(Stencil.D3Q19), LBStencil(Stencil.D3Q27)]:
        target = Target.GPU if ctx.cuda else Target.CPU
        data_type = "float64" if ctx.double_accuracy else "float32"
        openmp = True if ctx.openmp else False
        if ctx.optimize_for_localhost:
            cpu_vec = {"nontemporal": False, "assume_aligned": True}
        else:
            cpu_vec = None

        method = Method.SRT
        relaxation_rate = sp.symbols("omega")
        streaming_pattern = 'pull'

        pdfs = fields(f"pdfs({stencil.Q}): {data_type}[{stencil.D}D]", layout='fzyx')
        density_field, velocity_field = fields(f"density(1), velocity({stencil.D}): {data_type}[{stencil.D}D]",
                                               layout='fzyx')

        macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

        lbm_config = LBMConfig(stencil=stencil, method=method, relaxation_rate=relaxation_rate,
                               streaming_pattern=streaming_pattern)
        lbm_opt = LBMOptimisation(cse_global=False)

        collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

        generate_lbm_sweep_collection(ctx, f'{stencil.name}{method.name}', collision_rule,
                                      streaming_pattern='pull',
                                      field_layout='zyxf',
                                      refinement_scaling=None,
                                      macroscopic_fields=macroscopic_fields,
                                      target=target, data_type=data_type,
                                      cpu_openmp=openmp, cpu_vectorize_info=cpu_vec)

        ctx.write_all_files()
