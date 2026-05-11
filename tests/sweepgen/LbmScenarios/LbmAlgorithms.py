import argparse

import sympy as sp
from pystencilssfg import SourceFileGenerator

from pystencils import fields, Field, Target
from lbmpy import (
    LBStencil,
    Stencil,
    LBMConfig,
    LBMOptimisation,
    Method,
    ForceModel,
)

from sweepgen import Sweep, get_build_config
from sweepgen.prefabs import LbmBulk
from sweepgen.boundaries.grid_aligned import GridAlignedFreeSlip, GridAlignedNoSlip
from sweepgen.boundaries import SparseBoundary, FreeSlip

from sweepgen.build_config import DEBUG

DEBUG.use_hip_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    Sweep.use_v8core_fields()
    sfg.namespace("gen")

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--language", dest="language", default="CXX", type=str)

    args = parser.parse_args(sfg.context.argv)

    build_config = get_build_config(sfg)

    match args.language:
        case "CXX":
            target = Target.CurrentCPU
        case "HIP":
            target = Target.HIP
        case "CUDA":
            target = Target.CUDA
        case _:
            raise ValueError(f"Unexpected target id: {args.target}")

    build_config.override.target = target

    stencil = LBStencil(Stencil.D3Q19)
    d, q = stencil.D, stencil.Q
    f: Field
    f_tmp: Field
    f, f_tmp, rho, u = fields(
        f"f({q}), f_tmp({q}), rho, u({d}): double[{d}D]", layout="fzyx"
    )  # type: ignore
    force = sp.symbols(f"F_:{d}")

    stencil_name = stencil.name
    sfg.include(f"stencil/{stencil_name}.h")
    sfg.code(f"using LbStencil = walberla::stencil::{stencil_name};")

    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.TRT,
        output={"density": rho, "velocity": u},
        force=force,
        force_model=ForceModel.HE,
    )
    lbm_opt = LBMOptimisation(
        symbolic_field=f,
        symbolic_temporary_field=f_tmp,
    )

    lb_bulk = LbmBulk(
        sfg,
        "LBM",
        lbm_config,
        lbm_opt
    )

    sfg.generate(lb_bulk)

    with sfg.namespace("bc_grid_aligned"):
        sfg.generate(GridAlignedNoSlip("NoSlipTop", lb_bulk.lb_method, f, wall_orientation=(0, 0, 1)))

        sfg.generate(GridAlignedNoSlip("NoSlipBottom", lb_bulk.lb_method, f, wall_orientation=(0, 0, -1)))

        sfg.generate(
            GridAlignedFreeSlip("FreeSlipTop", lb_bulk.lb_method, f, wall_orientation=(0, 0, 1))
        )

        sfg.generate(
            GridAlignedFreeSlip("FreeSlipBottom", lb_bulk.lb_method, f, wall_orientation=(0, 0, -1))
        )

    with sfg.namespace("bc_sparse"):
        sfg.generate(
            SparseBoundary(
                FreeSlip(),
                lb_bulk.lb_method,
                lb_bulk.pdfs
            )
        )
