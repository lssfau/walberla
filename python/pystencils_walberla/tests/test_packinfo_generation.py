import unittest

import pystencils as ps
from pystencils_walberla import generate_pack_info_for_field, generate_pack_info_from_kernel
from pystencils_walberla.cmake_integration import ManualCodeGenerationContext


class PackinfoGenTest(unittest.TestCase):

    @staticmethod
    def test_packinfo_walberla_gen():
        for openmp in (False, True):
            for da in (False, True):
                with ManualCodeGenerationContext(openmp=openmp, double_accuracy=da) as ctx:
                    dtype = "float64" if ctx.double_accuracy else "float32"
                    f, g = ps.fields(f"f, g(4): {dtype}[3D]")
                    generate_pack_info_for_field(ctx, 'PI1', f)

                    src, dst = ps.fields(f"src, src_tmp: {dtype}[2D]")
                    stencil = [[0, -1, 0],
                               [-1, 4, -1],
                               [0, -1, 0]]
                    assignments = [ps.assignment_from_stencil(stencil, src, dst, normalization_factor=4)]
                    generate_pack_info_from_kernel(ctx, 'PI2', assignments)

                    expected_files = ('PI1.cpp', 'PI1.h', 'PI2.cpp', 'PI2.h')
                    assert all(e in ctx.files for e in expected_files)

                    for file_name_to_test in ('PI1.cpp', 'PI2.cpp'):
                        file_to_test = ctx.files[file_name_to_test]

                        # For Packing kernels it is better to not have OpenMP in the code because the packing kernels
                        # themselves are launched in an OpenMP environment. Still it could be forced but setting
                        # openmp to True in generate_pack_info_from_kernel
                        if openmp:
                            assert '#pragma omp parallel' not in file_to_test

                        if da:
                            assert 'float ' not in file_to_test
                        else:
                            assert 'double ' not in file_to_test


if __name__ == '__main__':
    unittest.main()
