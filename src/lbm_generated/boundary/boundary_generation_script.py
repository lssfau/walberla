import sympy as sp

from pystencils import Target

from lbmpy.creationfunctions import create_lb_method
from lbmpy import LBMConfig, Stencil, Method, LBStencil
from lbmpy.boundaries import ExtrapolationOutflow, FixedDensity, FreeSlip, NoSlip, UBB

from pystencils_walberla import ManualCodeGenerationContext, generate_info_header
from lbmpy_walberla.boundary_collection import generate_boundary_collection
from lbmpy_walberla import lbm_boundary_generator

with ManualCodeGenerationContext(openmp=False, optimize_for_localhost=False,
                                 mpi=True, double_accuracy=True, cuda=False) as ctx:

    for stencil in [LBStencil(Stencil.D3Q19), LBStencil(Stencil.D3Q27)]:
        target = Target.GPU if ctx.cuda else Target.CPU
        data_type = "float64" if ctx.double_accuracy else "float32"

        method = Method.SRT
        relaxation_rate = sp.symbols("omega")
        streaming_pattern = 'pull'

        lbm_config = LBMConfig(stencil=stencil, method=method, relaxation_rate=relaxation_rate,
                               streaming_pattern=streaming_pattern)

        lb_method = create_lb_method(lbm_config=lbm_config)

        outflow_west_boundary = ExtrapolationOutflow(normal_direction=(1, 0, 0), lb_method=lb_method)
        fixed_density_boundary = FixedDensity(density=sp.Symbol("density"))
        free_slip_boundary = FreeSlip(stencil)
        no_slip_boundary = NoSlip()
        ubb_boundary = UBB(sp.symbols("u_x, u_y, u_z"), data_type=data_type)

        outflow = lbm_boundary_generator(class_name=f'Outflow{stencil.name}', flag_uid='Outflow',
                                         boundary_object=outflow_west_boundary)

        fixed_density = lbm_boundary_generator(class_name=f'FixedDensity{stencil.name}', flag_uid='FixedDensity',
                                               boundary_object=fixed_density_boundary)

        free_slip = lbm_boundary_generator(class_name=f'FreeSlip{stencil.name}', flag_uid='FreeSlip',
                                           boundary_object=free_slip_boundary)

        no_slip = lbm_boundary_generator(class_name=f'NoSlip{stencil.name}', flag_uid='NoSlip',
                                         boundary_object=no_slip_boundary)

        ubb = lbm_boundary_generator(class_name=f'UBB{stencil.name}', flag_uid='UBB',
                                     boundary_object=ubb_boundary)

        boundaries = [outflow, fixed_density, free_slip, no_slip, ubb]
        generate_boundary_collection(ctx, f'{stencil.name}BoundaryCollection', boundary_generators=boundaries,
                                     lb_method=lb_method, streaming_pattern=streaming_pattern,
                                     target=target)

        ctx.write_all_files()
