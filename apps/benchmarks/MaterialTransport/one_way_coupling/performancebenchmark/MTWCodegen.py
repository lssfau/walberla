import copy
import sympy as sp
import pystencils as ps
from sympy.core.add import Add
from sympy.codegen.ast import Assignment
import sys

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.partially_saturated_cells import PSMConfig

from lbmpy.boundaries import NoSlip, UBB, FixedDensity, FreeSlip, DiffusionDirichlet, NeumannByCopy
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
from lbmpy_walberla import generate_boundary

info_header = """
const char * infoStencil_fluid = "{stencil}";
const char * infoStencil_concentration = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil_fluid = LBStencil(Stencil.D3Q27)
    stencil_concentration = LBStencil(Stencil.D3Q27)
    omega = sp.Symbol("omega")  # for now same for both the sweeps
    init_density_fluid = sp.Symbol("init_density_fluid")
    init_density_concentration = sp.Symbol("init_density_concentration")
    init_velocity_fluid = sp.symbols("init_velocity_fluid_:3")
    #init_velocity_concentration = sp.symbols("init_velocity_concentration_:3")
    pdfs_inter_fluid = sp.symbols("pdfs_inter_fluid:" + str(stencil_fluid.Q))
    pdfs_inter_concentration = sp.symbols("pdfs_inter_concentration:" + str(stencil_concentration.Q))
    layout = "fzyx"
    config_tokens = ctx.config.split("_")
    print(config_tokens[0]," ", config_tokens[1])
    MaxParticlesPerCell = int(2)
    methods = {
        "srt": Method.SRT,
        "trt": Method.TRT,
        "mrt": Method.MRT,
        "cumulant": Method.MONOMIAL_CUMULANT,
        "srt-smagorinsky": Method.SRT,
        "trt-smagorinsky": Method.TRT,
    }

# Fluid PDFs and fields
    pdfs_fluid, pdfs_fluid_tmp, velocity_field, density_field = ps.fields(
        f"pdfs_fluid({stencil_fluid.Q}), pdfs_fluid_tmp({stencil_fluid.Q}), velocity_field({stencil_fluid.D}), density_field({1}): {data_type}[3D]",
        layout=layout,
    )

    # Concentration PDFs and fields
    pdfs_concentration, pdfs_concentration_tmp, concentration_field = ps.fields(
        f"pdfs_concentration({stencil_concentration.Q}), pdfs_concentration_tmp({stencil_concentration.Q}), concentration_field({1}): {data_type}[3D]",
        layout=layout,
    )

    # Determine the output based on the coupling mode

    if config_tokens[1]== "1":
        concentration_output = None
        force_on_fluid = sp.symbols("F_:3")
        print("One-way fluid-concentration coupling set")

    elif config_tokens[1] == "2":
        concentration_output = {"density": concentration_field}
        force_on_fluid = sp.Matrix([0, 0, 0])
        print("Two-way fluid-concentration coupling set")

    # Ensure force_on_fluid is defined in all paths before using it

    # Fluid LBM optimisation
    lbm_fluid_opt = LBMOptimisation(
        cse_global=True,
        symbolic_field=pdfs_fluid,
        symbolic_temporary_field=pdfs_fluid_tmp,
        field_layout=layout,
    )

    # Concentration LBM optimisation
    lbm_concentration_opt = LBMOptimisation(
        cse_global=True,
        symbolic_field=pdfs_concentration,
        symbolic_temporary_field=pdfs_concentration_tmp,
        field_layout=layout,
    )

    # Fluid LBM config
    lbm_fluid_config = LBMConfig(
        stencil=stencil_fluid,
        method=Method.SRT,
        relaxation_rate=omega,
        output={"velocity": velocity_field},
        force= force_on_fluid,
        force_model=ForceModel.LUO,
        compressible=True,
    )

    # Concentration LBM config
    lbm_concentration_config = LBMConfig(
        stencil=stencil_concentration,
        method=Method.SRT,
        relaxation_rate=omega,
        velocity_input=velocity_field,
        #output={"density": concentration_output},
        force=sp.symbols("F_:3"),
        force_model=ForceModel.LUO,
        compressible=True,
        zero_centered=False
    )

    # =====================
    # Generate method
    # =====================

    method_fluid = create_lb_method(lbm_config=lbm_fluid_config)
    method_concentration = create_lb_method(lbm_config=lbm_concentration_config)

    pdfs_fluid_setter = macroscopic_values_setter(
        method_fluid, density=init_density_fluid, velocity=velocity_field.center_vector, pdfs=pdfs_fluid.center_vector
    )

    pdfs_concentration_setter = macroscopic_values_setter(
        method_concentration, density=concentration_field.center, velocity= velocity_field.center_vector,pdfs=pdfs_concentration.center_vector
    )

    # specify the target

    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU

    # Generate files

    generate_sweep(
        ctx,
        "LBMFluidSweep",
        create_lb_update_rule(lbm_config=lbm_fluid_config, lbm_optimisation=lbm_fluid_opt),
        field_swaps=[(pdfs_fluid, pdfs_fluid_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "LBMConcentrationSweep",
        create_lb_update_rule(lbm_config=lbm_concentration_config, lbm_optimisation=lbm_concentration_opt),
        field_swaps=[(pdfs_concentration, pdfs_concentration_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "LBMFluidSplitSweep",
        create_lb_update_rule(lbm_config=lbm_fluid_config, lbm_optimisation=lbm_fluid_opt),
        field_swaps=[(pdfs_fluid, pdfs_fluid_tmp)],
        target=target,
        inner_outer_split=True,
    )

    generate_sweep(
        ctx,
        "LBMConcentrationSplitSweep",
        create_lb_update_rule(lbm_config=lbm_concentration_config, lbm_optimisation=lbm_concentration_opt),
        field_swaps=[(pdfs_concentration, pdfs_concentration_tmp)],
        target=target,
        inner_outer_split=True,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PackInfoFluid",
        create_lb_update_rule(lbm_config=lbm_fluid_config, lbm_optimisation=lbm_fluid_opt),
        target=target,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PackInfoConcentration",
        create_lb_update_rule(lbm_config=lbm_concentration_config, lbm_optimisation=lbm_concentration_opt),
        target=target,
    )

    generate_sweep(ctx, "InitializeFluidDomain", pdfs_fluid_setter, target=target)
    generate_sweep(ctx, "InitializeConcentrationDomain", pdfs_concentration_setter, target=target)

    # Fluid Boundary conditions
    generate_boundary(
        ctx,
        "BC_Fluid_NoSlip",
        NoSlip(),
        method_fluid,
        field_name=pdfs_fluid.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_velocity_fluid = sp.symbols("bc_velocity_fluid_:3")
    generate_boundary(
        ctx,
        "BC_Fluid_UBB",
        UBB(bc_velocity_fluid),
        method_fluid,
        field_name=pdfs_fluid.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_density_fluid = sp.Symbol("bc_density_fluid")
    generate_boundary(
        ctx,
        "BC_Fluid_Density",
        FixedDensity(bc_density_fluid),
        method_fluid,
        field_name=pdfs_fluid.name,
        streaming_pattern="pull",
        target=target,
    )

    generate_boundary(
        ctx,
        "BC_Fluid_FreeSlip",
        FreeSlip(stencil_fluid),
        method_fluid,
        field_name=pdfs_fluid.name,
        streaming_pattern="pull",
        target=target,
    )


    # Concentration Boundary conditions

    generate_boundary(
        ctx,
        "BC_Concentration_NoSlip",
        NoSlip(),
        method_concentration,
        field_name=pdfs_concentration.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_velocity_concentration = sp.symbols("bc_velocity_concentration_:3")   ## is it needed ?
    generate_boundary(
        ctx,
        "BC_Concentration_UBB",
        UBB(bc_velocity_concentration),
        method_concentration,
        field_name=pdfs_concentration.name,
        streaming_pattern="pull",
        target=target,
    )

    bc_density_concentration = sp.Symbol("bc_density_concentration")
    generate_boundary(
        ctx,
        "BC_Concentration_Density",
        DiffusionDirichlet(bc_density_concentration),
        method_concentration,
        field_name=pdfs_concentration.name,
        streaming_pattern="pull",
        target=target,
    )

    generate_boundary(
        ctx,
        "BC_Concentration_FreeSlip",
        FreeSlip(stencil_concentration),
        method_concentration,
        field_name=pdfs_concentration.name,
        streaming_pattern="pull",
        target=target,
    )

    generate_boundary(
        ctx,
        "BC_Concentration_Neumann",
        NeumannByCopy(stencil_concentration),
        method_concentration,
        field_name=pdfs_concentration.name,
        streaming_pattern="pull",
        target=target,
    )



    # Info header containing correct template definitions for stencil and fields
    infoHeaderParams = {
        "stencil_fluid": stencil_fluid.name,
        "stencil_concentration": stencil_concentration.name,
        "streaming_pattern": lbm_fluid_config.streaming_pattern,
        "collision_setup": Method.SRT,
        "cse_global": int(lbm_fluid_opt.cse_global),
        "cse_pdfs": int(lbm_fluid_opt.cse_pdfs),
    }

    additional_code = f"""
    const char * infoStencil_concentration = "{stencil_concentration.name}";
    const char * infoStencil_fluid = "{stencil_fluid.name}";
    const char * infoStreamingPattern = "{lbm_fluid_config.streaming_pattern}";
    const char * infoCollisionSetup = "{Method.SRT}";
    const bool infoCseGlobal = {int(lbm_fluid_opt.cse_global)};
    const bool infoCsePdfs = {int(lbm_fluid_opt.cse_pdfs)};
    """


    stencil_typedefs = {"Stencil_Fluid_T": stencil_fluid, "CommunicationStencil_Fluid_T": stencil_fluid, "Stencil_Concentration_T": stencil_concentration, "CommunicationStencil_Concentration_T": stencil_concentration}
    field_typedefs = {
        "PdfField_fluid_T": pdfs_fluid,
        "DensityField_fluid_T": density_field,
        "VelocityField_fluid_T": velocity_field,
        "PdfField_concentration_T": pdfs_concentration,
        "DensityField_concentration_T": concentration_field,
    }

    generate_info_header(
        ctx,
        "GeneralInfoHeader",
        stencil_typedefs=stencil_typedefs,
        field_typedefs=field_typedefs,
        #additional_code=additional_code,
    )

    # Getter & setter to compute moments from pdfs
    pdfs_fluid_getter = macroscopic_values_getter(
        method_fluid, density=density_field, velocity=velocity_field.center_vector,pdfs=pdfs_fluid.center_vector
    )

    pdfs_concentration_getter = macroscopic_values_getter(
        method_concentration, density=concentration_field, velocity=None,pdfs=pdfs_concentration.center_vector
    )

    generate_sweep(ctx, "FluidMacroSetter", pdfs_fluid_setter)
    generate_sweep(ctx, "FluidMacroGetter", pdfs_fluid_getter)
    generate_sweep(ctx, "ConcentrationMacroSetter", pdfs_concentration_setter)
    generate_sweep(ctx, "ConcentrationMacroGetter", pdfs_concentration_getter)