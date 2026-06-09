import pystencils as ps

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.build_config import DEBUG
from sweepgen.symbolic import cell, cell_index

DEBUG.use_cpu_default()

with SourceFileGenerator() as sfg:
    get_build_config(sfg).override.target = ps.Target.GenericCPU

    out_real = ps.fields("out(3): double[3D]")
    out_int = ps.fields("out(3): int64_t[3D]")

    #   Global Cell Centers

    asms = [
        ps.Assignment(out_real.center(0), cell.x()),
        ps.Assignment(out_real.center(1), cell.y()),
        ps.Assignment(out_real.center(2), cell.z()),
    ]

    sweep = Sweep("CellCentersGlobal", asms)
    sfg.generate(sweep)

    #   Local Cell Centers

    @ps.flow.block
    def cellCentersLocal(_eq):
        _eq.store[out_real.center(0)] = cell.local_x()
        _eq.store[out_real.center(1)] = cell.local_y()
        _eq.store[out_real.center(2)] = cell.local_z()

    sweep = Sweep("CellCentersLocal", cellCentersLocal)
    sfg.generate(sweep)

    #   Local Cell Indices

    @ps.flow.block
    def cellIndicesLocal(_eq):
        _eq.store[out_int.center(0)] = cell_index.x_local()
        _eq.store[out_int.center(1)] = cell_index.y_local()
        _eq.store[out_int.center(2)] = cell_index.z_local()

    sweep = Sweep("CellIndicesLocal", cellIndicesLocal)
    sfg.generate(sweep)

    #   Global Cell Indices

    asms = [
        ps.Assignment(out_int.center(0), cell_index.x_global()),
        ps.Assignment(out_int.center(1), cell_index.y_global()),
        ps.Assignment(out_int.center(2), cell_index.z_global()),
    ]

    sweep = Sweep("CellIndicesGlobal", asms)
    sfg.generate(sweep)
