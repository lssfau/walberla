import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"
    # ----- Solving the 2D Poisson equation with rhs --------------------------
    dx = sp.Symbol("dx")
    dy = sp.Symbol("dy")
    src, dst, rhs = ps.fields(f"src, src_tmp, rhs: {field_type}[2D]", layout='fzyx')

    @ps.kernel
    def kernel_func():
        dst[0, 0] @= ((dy**2 * (src[1, 0] + src[-1, 0]))
                      + (dx**2 * (src[0, 1] + src[0, -1]))
                      - (rhs[0, 0] * dx**2 * dy**2)) / (2 * (dx**2 + dy**2))

    generate_sweep(ctx, 'PoissonGPU', kernel_func, field_swaps=[(src, dst), ], target=ps.Target.GPU)
