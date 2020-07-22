import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.boundaries import NoSlip, UBB
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import RefinementScaling, generate_boundary, generate_lattice_model

with CodeGeneration() as ctx:
    omega, omega_free = sp.symbols("omega, omega_free")
    force_field, vel_field, omega_out = ps.fields("force(3), velocity(3), omega_out: [3D]", layout='zyxf')

    # the collision rule of the LB method where the some advanced features
    collision_rule = create_lb_collision_rule(
        stencil='D3Q19', compressible=True,
        method='mrt', relaxation_rates=[omega, omega, omega_free, omega_free, omega_free, omega_free],
        entropic=True,                    # entropic method where second omega is chosen s.t. entropy condition
        omega_output_field=omega_out,     # scalar field where automatically chosen omega of entropic or
                                          # Smagorinsky method is written to
        force=force_field.center_vector,  # read forces for each lattice cell from an external force field
                                          # that is initialized and changed in C++ app
        output={'velocity': vel_field},   # write macroscopic velocity to field in every time step
                                          # useful for coupling multiple LB methods,
                                          # e.g. hydrodynamic to advection/diffusion LBM
        optimization={'cse_global': True}
    )

    # the refinement scaling object describes how certain parameters are scaled across grid scales
    # there are two default scaling behaviors available for relaxation rates and forces:
    scaling = RefinementScaling()
    scaling.add_standard_relaxation_rate_scaling(omega)
    scaling.add_force_scaling(force_field)

    # generate lattice model and (optionally) boundary conditions
    # for CPU simulations waLBerla's internal boundary handling can be used as well
    generate_lattice_model(ctx, 'LbCodeGenerationExample_LatticeModel', collision_rule, refinement_scaling=scaling)
    generate_boundary(ctx, 'LbCodeGenerationExample_UBB', UBB([0.05, 0, 0]), collision_rule.method)
    generate_boundary(ctx, 'LbCodeGenerationExample_NoSlip', NoSlip(), collision_rule.method)
