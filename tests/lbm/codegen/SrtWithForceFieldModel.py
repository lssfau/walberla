import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method
from lbmpy.boundaries import NoSlip, UBB
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model, RefinementScaling, generate_boundary


with CodeGeneration() as ctx:
    omega = sp.Symbol("omega")
    force_field = ps.fields("force(3): [3D]", layout='fzyx')

    # lattice Boltzmann method
    lb_method = create_lb_method(stencil='D3Q19', method='srt', relaxation_rates=[omega],
                                 force_model='guo', force=force_field.center_vector)

    scaling = RefinementScaling()
    scaling.add_standard_relaxation_rate_scaling(omega)
    scaling.add_force_scaling(force_field)

    # generate components
    generate_lattice_model(ctx, 'SrtWithForceFieldModel', lb_method, refinement_scaling=scaling)
    generate_boundary(ctx, 'MyUBB', UBB([0.05, 0, 0]), lb_method)
    generate_boundary(ctx, 'MyNoSlip', NoSlip(), lb_method)
