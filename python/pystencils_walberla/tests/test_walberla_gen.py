import unittest

import sympy as sp

import pystencils as ps
from pystencils_walberla import generate_sweep
from pystencils_walberla.cmake_integration import ManualCodeGenerationContext


class CodegenTest(unittest.TestCase):

    @staticmethod
    def test_codegen():
        for openmp in (False, True):
            for da in (False, True):
                with ManualCodeGenerationContext(openmp=openmp, double_accuracy=da) as ctx:
                    h = sp.symbols("h")

                    dtype = "float64" if ctx.double_accuracy else "float32"

                    # ----- Jacobi 2D - created by specifying weights in nested list --------------------------
                    src, dst = ps.fields(f"src, src_tmp: {dtype}[2D]")
                    stencil = [[0, -1, 0],
                               [-1, 4, -1],
                               [0, -1, 0]]
                    assignments = ps.assignment_from_stencil(stencil, src, dst, normalization_factor=4 * h ** 2)
                    generate_sweep(ctx, 'JacobiKernel2D', assignments, field_swaps=[(src, dst)])

                    # ----- Jacobi 3D - created by using kernel_decorator with assignments in '@=' format -----
                    src, dst = ps.fields(f"src, src_tmp: {dtype}[3D]")

                    @ps.kernel
                    def kernel_func():
                        dst[0, 0, 0] @= (src[1, 0, 0] + src[-1, 0, 0]
                                         + src[0, 1, 0] + src[0, -1, 0]
                                         + src[0, 0, 1] + src[0, 0, -1]) / (6 * h ** 2)

                    generate_sweep(ctx, 'JacobiKernel3D', kernel_func, field_swaps=[(src, dst)])

                    expected_files = ('JacobiKernel3D.cpp', 'JacobiKernel3D.h',
                                      'JacobiKernel2D.cpp', 'JacobiKernel2D.h')
                    assert all(e in ctx.files for e in expected_files)

                    for file_name_to_test in ('JacobiKernel3D.cpp', 'JacobiKernel2D.cpp'):
                        file_to_test = ctx.files[file_name_to_test]
                        if openmp:
                            assert '#pragma omp parallel' in file_to_test

                        if da:
                            assert 'float ' not in file_to_test
                        else:
                            assert 'double ' not in file_to_test


if __name__ == '__main__':
    unittest.main()
