import pystencils as ps
from pystencilssfg import SourceFileGenerator, SfgComposer

from sweepgen import Sweep
from sweepgen.symbolic import cell
from sweepgen.build_config import DEBUG, get_build_config

DEBUG.use_cpu_default()


def test_1d_2d_fields(sfg: SfgComposer):
    coords = [cell.x(), cell.y(), cell.z()]

    for d in (1, 2):
        with sfg.namespace(f"Kernels{d}D"):
            for layout in ("fzyx", "zyxf"):
                rhs = coords
                f = ps.fields(f"f(3): double[{d}D]", layout=layout)

                asms = [ps.Assignment(f(i), r) for i, r in enumerate(rhs)]
                sfg.generate(Sweep(f"Set_{layout}", asms))


def test_experimental_fields(sfg: SfgComposer):
    with sfg.namespace("TestExperimentalFields"):
        with Sweep.use_v8core_fields():
            f, g = ps.fields("f(3), g: double[3D]")

            init_asms = [
                ps.Assignment(f(i), c)
                for i, c in enumerate([cell.x(), cell.y(), cell.z()])
            ]
            sfg.generate(Sweep("InitializeF", init_asms))

            sum_asm = ps.Assignment(g(), f(0) + f(1) + f(2))
            sfg.generate(Sweep("PointwiseSum", [sum_asm]))


def test_1d_2d_fields_v8(sfg: SfgComposer):
    coords = [cell.x(), cell.y(), cell.z()]

    with Sweep.use_v8core_fields():
        for d in (1, 2):
            with sfg.namespace(f"v8::Kernels{d}D"):
                for layout in ("fzyx", "zyxf"):
                    rhs = coords
                    f = ps.fields(f"f(3): double[{d}D]", layout=layout)

                    asms = [ps.Assignment(f(i), r) for i, r in enumerate(rhs)]
                    sfg.generate(Sweep(f"Set_{layout}", asms))


def test_field_swaps(sfg: SfgComposer):
    with sfg.namespace("TestFieldSwaps"):
        with Sweep.use_v8core_fields():
            f, f_tmp = ps.fields("f, f_tmp: double[3D]")

            stencil_asm = ps.Assignment(
                f_tmp(),
                f[-1, 0, 0]
                + f[1, 0, 0]
                + f[0, -1, 0]
                + f[0, 1, 0]
                + f[0, 0, -1]
                + f[0, 0, 1],
            )

            sweep = Sweep("TestStencil", [stencil_asm])
            sweep.swap_fields(f, f_tmp)
            sfg.generate(sweep)


with SourceFileGenerator() as sfg:
    get_build_config(sfg).target = ps.Target.CPU

    with sfg.namespace("gen"):
        test_1d_2d_fields(sfg)
        test_experimental_fields(sfg)
        test_1d_2d_fields_v8(sfg)
        test_field_swaps(sfg)
