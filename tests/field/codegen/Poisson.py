import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    # ----- Solving the 2D Poisson equation with constant rhs --------------------------
    rhs = sp.Symbol("rhs")
    h = sp.Symbol("h")
    src, dst = ps.fields("src, src_tmp: [2D]", layout='fzyx')

    @ps.kernel
    def kernel_func():
        dst[0, 0] @= (src[1, 0] + src[-1, 0] +
                      src[0, 1] + src[0, -1]) / 4 + ((rhs * h**2) / 4)
    generate_sweep(ctx, 'Poisson', kernel_func, field_swaps=[(src, dst)])
