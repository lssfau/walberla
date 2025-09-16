import pystencils as ps
from pystencilssfg import SourceFileGenerator

from sweepgen import Sweep
from sweepgen.symbolic import cell
from sweepgen.build_config import DEBUG, get_build_config

DEBUG.use_cpu_default()

with SourceFileGenerator() as sfg:
    get_build_config(sfg).target = ps.Target.CPU

    coords = [cell.x(), cell.y(), cell.z()]

    for d in (1, 2):
        with sfg.namespace(f"Kernels{d}D"):
            for layout in ("fzyx", "zyxf"):
                rhs = coords
                f = ps.fields(f"f(3): double[{d}D]", layout=layout)

                asms = [ps.Assignment(f(i), r) for i, r in enumerate(rhs)]
                sfg.generate(Sweep(f"Set_{layout}", asms))
