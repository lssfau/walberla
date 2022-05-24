import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


def grad(f):
    return sp.Matrix([ps.fd.diff(f, i) for i in range(f.spatial_dimensions)])


D = sp.Symbol("D")
        
with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"

    c = ps.fields(f"c: {data_type}[2D]", layout='fzyx')
    j = ps.fields(f"j(2): {data_type}[2D]", layout='fzyx', field_type=ps.FieldType.STAGGERED_FLUX)

    ek = ps.fd.FVM1stOrder(c, flux=-D * grad(c))

    generate_sweep(ctx, 'EKFlux', ek.discrete_flux(j), staggered=True)
    generate_sweep(ctx, 'EKContinuity', ek.discrete_continuity(j))
