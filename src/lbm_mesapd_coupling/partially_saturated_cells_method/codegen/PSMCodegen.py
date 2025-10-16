import copy
import sympy as sp
import pystencils as ps
from sympy.core.add import Add
from sympy.codegen.ast import Assignment

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.partially_saturated_cells import PSMConfig

from lbmpy.boundaries import NoSlip, UBB, FixedDensity, FreeSlip
from lbmpy.creationfunctions import (
    create_lb_update_rule,
    create_lb_method,
    create_psm_update_rule,
)

from lbmpy.macroscopic_value_kernels import (
    macroscopic_values_getter,
    macroscopic_values_setter,
)

from pystencils_walberla import (
    CodeGeneration,
    generate_info_header,
    generate_sweep,
    generate_pack_info_from_kernel,
)

from pystencils_walberla.utility import get_vectorize_instruction_set

from lbmpy_walberla import generate_boundary

# Based on the following paper: https://doi.org/10.1016/j.compfluid.2017.05.033

info_header = """
const char * infoStencil = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    if ctx.optimize_for_localhost:
        isa = get_vectorize_instruction_set(ctx)
        if isa in ("neon", "sve", "sve2", "sme"):
            ctx.optimize_for_localhost = False

    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)
    omega = sp.Symbol("omega")
    init_density = sp.Symbol("init_density")
    init_velocity = sp.symbols("init_velocity_:3")
    pdfs_inter = sp.symbols("pdfs_inter:" + str(stencil.Q))
    layout = "fzyx"
    config_tokens = ctx.config.split("_")
    MaxParticlesPerCell = int(config_tokens[2])
    methods = {
        "srt": Method.SRT,
        "trt": Method.TRT,
        "mrt": Method.MRT,
        "cumulant": Method.MONOMIAL_CUMULANT,
        "srt-smagorinsky": Method.SRT,
        "trt-smagorinsky": Method.TRT,
    }
    # Solid collision variant
    SC = int(config_tokens[1][2])

    pdfs, pdfs_tmp, velocity_field, density_field = ps.fields(
        f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}), velocity_field({stencil.D}), density_field({1}): {data_type}[3D]",
        layout=layout,
    )

    particle_velocities, particle_forces, Bs = ps.fields(
        f"particle_v({MaxParticlesPerCell * stencil.D}), particle_f({MaxParticlesPerCell * stencil.D}), Bs({MaxParticlesPerCell}): {data_type}[3D]",
        layout=layout,
    )

    # Solid fraction field
    B = ps.fields(f"b({1}): {data_type}[3D]", layout=layout)

    psm_opt = LBMOptimisation(
        cse_global=True,
        symbolic_field=pdfs,
        symbolic_temporary_field=pdfs_tmp,
        field_layout=layout,
    )

    psm_config = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        SC=SC,
        MaxParticlesPerCell=MaxParticlesPerCell,
        individual_fraction_field=Bs,
        particle_force_field=particle_forces,
    )

    lbm_config = LBMConfig(
        stencil=stencil,
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        force=sp.symbols("F_:3"),
        force_model=ForceModel.LUO,
        compressible=True,
        psm_config=psm_config,
    )

    if config_tokens[0] == "srt-smagorinsky" or config_tokens[0] == "trt-smagorinsky":
        lbm_config.smagorinsky = True

    # =====================
    # Generate method
    # =====================

    method = create_lb_method(lbm_config=lbm_config)

    node_collection = create_psm_update_rule(lbm_config, psm_opt)

    pdfs_setter = macroscopic_values_setter(
        method, init_density, init_velocity, pdfs.center_vector
    )

    # Use average velocity of all intersecting particles when setting PDFs (mandatory for SC=3)
    for i, sub_exp in enumerate(pdfs_setter.subexpressions[-3:]):
        rhs = []
        for summand in sub_exp.rhs.args:
            rhs.append(summand * (1.0 - B.center))
        for p in range(MaxParticlesPerCell):
            rhs.append(particle_velocities(p * stencil.D + i) * Bs.center(p))
        pdfs_setter.subexpressions.remove(sub_exp)
        pdfs_setter.subexpressions.append(Assignment(sub_exp.lhs, Add(*rhs)))

    # =====================
    # Write method to files
    # =====================

    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU

    # Generate files
    generate_sweep(
        ctx,
        "PSMSweep",
        node_collection,
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "PSMSweepSplit",
        node_collection,
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
        inner_outer_split=True,
    )

    config_without_psm = LBMConfig(
        stencil=stencil,
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        force=sp.symbols("F_:3"),
        force_model=ForceModel.LUO,
        compressible=True,
    )

    if config_tokens[0] == "srt-smagorinsky" or config_tokens[0] == "trt-smagorinsky":
        config_without_psm.smagorinsky = True

    generate_sweep(
        ctx,
        "LBMSweep",
        create_lb_update_rule(lbm_config=config_without_psm, lbm_optimisation=psm_opt),
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "LBMSplitSweep",
        create_lb_update_rule(lbm_config=config_without_psm, lbm_optimisation=psm_opt),
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
        inner_outer_split=True,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PSMPackInfo",
        create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=psm_opt),
        target=target,
    )

    generate_sweep(ctx, "InitializeDomainForPSM", pdfs_setter, target=target)

    # Boundary conditions
    generate_boundary(
        ctx,
        "PSM_NoSlip",
        NoSlip(),
        method,
        field_name=pdfs.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_velocity = sp.symbols("bc_velocity_:3")
    generate_boundary(
        ctx,
        "PSM_UBB",
        UBB(bc_velocity),
        method,
        field_name=pdfs.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_density = sp.Symbol("bc_density")
    generate_boundary(
        ctx,
        "PSM_Density",
        FixedDensity(bc_density),
        method,
        field_name=pdfs.name,
        streaming_pattern="pull",
        target=target,
    )

    generate_boundary(
        ctx,
        "PSM_FreeSlip",
        FreeSlip(stencil),
        method,
        field_name=pdfs.name,
        streaming_pattern="pull",
        target=target,
    )

    # Info header containing correct template definitions for stencil and fields
    infoHeaderParams = {
        "stencil": stencil.name,
        "streaming_pattern": lbm_config.streaming_pattern,
        "collision_setup": config_tokens[0],
        "cse_global": int(psm_opt.cse_global),
        "cse_pdfs": int(psm_opt.cse_pdfs),
    }

    stencil_typedefs = {"Stencil_T": stencil, "CommunicationStencil_T": stencil}
    field_typedefs = {
        "PdfField_T": pdfs,
        "DensityField_T": density_field,
        "VelocityField_T": velocity_field,
    }

    generate_info_header(
        ctx,
        "PSM_InfoHeader",
        stencil_typedefs=stencil_typedefs,
        field_typedefs=field_typedefs,
        additional_code=info_header.format(**infoHeaderParams),
    )

    # Getter & setter to compute moments from pdfs
    setter_assignments = macroscopic_values_setter(
        method,
        velocity=velocity_field.center_vector,
        pdfs=pdfs.center_vector,
        density=1.0,
    )
    getter_assignments = macroscopic_values_getter(
        method,
        density=density_field,
        velocity=velocity_field.center_vector,
        pdfs=pdfs.center_vector,
    )
    generate_sweep(ctx, "PSM_MacroSetter", setter_assignments)
    generate_sweep(ctx, "PSM_MacroGetter", getter_assignments)
