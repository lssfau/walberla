import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    h = sp.symbols("h")

    # ----- Jacobi 2D - created by specifying weights in nested list --------------------------
    src, dst = ps.fields("src, src_tmp: [2D]")
    stencil = [[0, -1, 0],
               [-1, 4, -1],
               [0, -1, 0]]
    assignments = ps.assignment_from_stencil(stencil, src, dst, normalization_factor=4 * h**2)
    generate_sweep(ctx, 'CudaJacobiKernel2D', assignments, field_swaps=[(src, dst)], target="gpu")

    # ----- Jacobi 3D - created by using kernel_decorator with assignments in '@=' format -----
    src, dst = ps.fields("src, src_tmp: [3D]")

    @ps.kernel
    def kernel_func():
        dst[0, 0, 0] @= (src[1, 0, 0] + src[-1, 0, 0] +
                         src[0, 1, 0] + src[0, -1, 0] +
                         src[0, 0, 1] + src[0, 0, -1]) / (6 * h ** 2)

    generate_sweep(ctx, 'CudaJacobiKernel3D', kernel_func, field_swaps=[(src, dst)], target="gpu")
