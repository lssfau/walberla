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
from lbmpy_walberla import generate_boundary

info_header = """
const char * infoStencil = "{stencil}";
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
    init_velocity_concentration = sp.symbols("init_velocity_concentration_:3")
    pdfs_inter_fluid = sp.symbols("pdfs_inter_fluid:" + str(stencil_fluid.Q))
    pdfs_inter_concentration = sp.symbols("pdfs_inter_concentration:" + str(stencil_concentration.Q))
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
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        force=sp.symbols("F_:3"),
        force_model=ForceModel.LUO,
        compressible=True,
    )

    # Concentration LBM config
    lbm_concentration_config = LBMConfig(
        stencil=stencil_concentration,
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        force=sp.symbols("F_:3"),
        force_model=ForceModel.LUO,
        compressible=True,
    )