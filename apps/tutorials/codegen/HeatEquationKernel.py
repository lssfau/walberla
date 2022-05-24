import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"

    u, u_tmp = ps.fields(f"u, u_tmp: {data_type}[2D]", layout='fzyx')
    kappa = sp.Symbol("kappa")
    dx = sp.Symbol("dx")
    dt = sp.Symbol("dt")
    heat_pde = ps.fd.transient(u) - kappa * (ps.fd.diff(u, 0, 0) + ps.fd.diff(u, 1, 1))

    discretize = ps.fd.Discretization2ndOrder(dx=dx, dt=dt)
    heat_pde_discretized = discretize(heat_pde)
    heat_pde_discretized = heat_pde_discretized.args[1] + heat_pde_discretized.args[0].simplify()


    @ps.kernel
    def update():
        u_tmp.center @= heat_pde_discretized


    ac = ps.AssignmentCollection(update)
    ac = ps.simp.simplifications.add_subexpressions_for_divisions(ac)

    generate_sweep(ctx, 'HeatEquationKernel', ac)
