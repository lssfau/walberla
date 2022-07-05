import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"
    # ----- Solving the 2D Poisson equation with rhs --------------------------
    dx2 = sp.Symbol("dx_square")
    dy2 = sp.Symbol("dy_square")
    src, dst, rhs = ps.fields(f"src, src_tmp, rhs: {field_type}[2D]", layout='fzyx')

    @ps.kernel
    def kernel_func():
        src[0, 0] @= ((dy2 * (src[1, 0] + src[-1, 0]))
                      + (dx2 * (src[0, 1] + src[0, -1]))
                      - (rhs[0, 0] * dx2 * dy2)) / (2.0 * (dx2 + dy2))

    generate_sweep(ctx, 'Poisson', kernel_func)
