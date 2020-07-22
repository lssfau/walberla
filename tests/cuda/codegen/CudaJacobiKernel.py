import numpy as np
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep


with CodeGeneration() as ctx:
    # ----- Stencil 2D - created by specifying weights in nested list --------------------------
    src, dst = ps.fields("src, src_tmp: [2D]", layout='fzyx')
    stencil = [[1.11, 2.22, 3.33],
               [4.44, 5.55, 6.66],
               [7.77, 8.88, 9.99]]
    assignments = ps.assignment_from_stencil(stencil, src, dst, normalization_factor=1 / np.sum(stencil))
    generate_sweep(ctx, 'CudaJacobiKernel2D', assignments, field_swaps=[(src, dst)], target="gpu")

    # ----- Stencil 3D - created by using kernel_decorator with assignments in '@=' format -----
    src, dst = ps.fields("src, src_tmp: [3D]")

    @ps.kernel
    def kernel_func():
        dst[0, 0, 0] @= (3 * src[1, 0, 0] + 4 * src[-1, 0, 0]
                         + 5 * src[0, 1, 0] + 6 * src[0, -1, 0]
                         + 7 * src[0, 0, 1] + 8 * src[0, 0, -1]) / 33

    generate_sweep(ctx, 'CudaJacobiKernel3D', kernel_func, field_swaps=[(src, dst)], target="gpu")
