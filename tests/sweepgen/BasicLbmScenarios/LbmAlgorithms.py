import argparse

import sympy as sp
from pystencilssfg import SourceFileGenerator

from pystencils import fields, Field, Target
from lbmpy import (
    LBStencil,
    Stencil,
    LBMConfig,
    LBMOptimisation,
    create_lb_update_rule,
    ForceModel,
)
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

from sweepgen import Sweep, get_build_config
from sweepgen.boundaries import NoSlip, FreeSlip
from sweepgen.communication import GpuPdfFieldPackInfo

from sweepgen.build_config import DEBUG

DEBUG.use_hip_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", dest="target", default="cpu", type=str)

    args = parser.parse_args(sfg.context.argv)

    build_config = get_build_config(sfg)

    match args.target:
        case "cpu":
            target = Target.CurrentCPU
            sfg.code("#define LBM_SCENARIOS_CPU_BUILD true")
        case "hip":
            target = Target.HIP
            sfg.code("#define LBM_SCENARIOS_GPU_BUILD true")
        case "cuda":
            target = Target.CUDA
            sfg.code("#define LBM_SCENARIOS_GPU_BUILD true")
        case _:
            raise ValueError(f"Unexpected target id: {args.target}")

    build_config.override.target = target

    sfg.namespace("BasicLbmScenarios::gen")

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
        output={"density": rho, "velocity": u},
        force=force,
        force_model=ForceModel.HE,
    )
    lbm_opt = LBMOptimisation(
        symbolic_field=f,
        symbolic_temporary_field=f_tmp,
    )

    lb_update = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    assert lb_update is not None

    with sfg.namespace("bulk"):
        lb_update_sweep = Sweep("LbStreamCollide", lb_update)
        lb_update_sweep.swap_fields(f, f_tmp)
        sfg.generate(lb_update_sweep)

        lb_init = macroscopic_values_setter(
            lb_update.method,
            density=rho.center,
            velocity=u.center_vector,
            pdfs=f,
            set_pre_collision_pdfs=True
        )
        lb_init_sweep = Sweep("LbInitFromFields", lb_init)
        sfg.generate(lb_init_sweep)

        lb_init_constant = macroscopic_values_setter(
            lb_update.method,
            density=sp.Symbol("density"),
            velocity=sp.symbols(f"velocity_:{d}"),
            pdfs=f,
            set_pre_collision_pdfs=True
        )
        lb_init_constant_sweep = Sweep("LbInitConstant", lb_init_constant)
        sfg.generate(lb_init_constant_sweep)

    with sfg.namespace("bc_grid_aligned"):
        sfg.generate(NoSlip("NoSlipTop", lb_update.method, f, wall_orientation=(0, 0, 1)))

        sfg.generate(NoSlip("NoSlipBottom", lb_update.method, f, wall_orientation=(0, 0, -1)))

        sfg.generate(
            FreeSlip("FreeSlipTop", lb_update.method, f, wall_orientation=(0, 0, 1))
        )

        sfg.generate(
            FreeSlip("FreeSlipBottom", lb_update.method, f, wall_orientation=(0, 0, -1))
        )

    with sfg.namespace("bc_sparse"):
        irreg_freeslip = FreeSlip(
            "FreeSlipIrregular",
            lb_update.method,
            f,
            wall_orientation=FreeSlip.IRREGULAR,
        )
        sfg.generate(irreg_freeslip)

    if build_config.override.target.is_gpu():
        with sfg.namespace("comm"):
            pack_info = GpuPdfFieldPackInfo("GpuPdfsPackInfo", stencil, f)
            sfg.generate(pack_info)
