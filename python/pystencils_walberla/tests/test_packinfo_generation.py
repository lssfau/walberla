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
                    f, g = ps.fields("f, g(4): {}[3D]".format(dtype))
                    generate_pack_info_for_field(ctx, 'PI1', f)

                    src, dst = ps.fields("src, src_tmp: {}[2D]".format(dtype))
                    stencil = [[0, -1, 0],
                               [-1, 4, -1],
                               [0, -1, 0]]
                    assignments = [ps.assignment_from_stencil(stencil, src, dst, normalization_factor=4)]
                    generate_pack_info_from_kernel(ctx, 'PI2', assignments)

                    expected_files = ('PI1.cpp', 'PI1.h', 'PI2.cpp', 'PI2.h')
                    assert all(e in ctx.files for e in expected_files)

                    for file_name_to_test in ('PI1.cpp', 'PI2.cpp'):
                        file_to_test = ctx.files[file_name_to_test]
                        if openmp:
                            assert '#pragma omp parallel' in file_to_test

                        if da:
                            assert 'float ' not in file_to_test
                        else:
                            assert 'double ' not in file_to_test


if __name__ == '__main__':
    unittest.main()
