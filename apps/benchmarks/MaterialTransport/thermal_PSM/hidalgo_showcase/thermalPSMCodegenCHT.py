import copy
import sympy as sp
import pystencils as ps
from sympy.core.add import Add
from sympy.codegen.ast import Assignment
import sys
import numpy as np

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.partially_saturated_cells import PSMConfig

from lbmpy.boundaries import NoSlip, UBB, FixedDensity, FreeSlip, DiffusionDirichlet, NeumannByCopy, SimpleExtrapolationOutflow,ExtrapolationOutflow
from lbmpy.creationfunctions import (
    create_lb_update_rule,
    create_lb_method,
    create_psm_update_rule,
    create_lb_collision_rule,
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
from lbmpy_walberla.additional_data_handler import DiffusionDirichletAdditionalDataHandler
from pystencils.cache import clear_cache
clear_cache()





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
    stencil_fluid = LBStencil(Stencil.D3Q19)
    stencil_concentration = LBStencil(Stencil.D3Q7)
    stencil_energy = LBStencil(Stencil.D3Q7)
    omega = sp.Symbol("omega")  # for now same for both the sweeps
    init_density_fluid = sp.Symbol("init_density_fluid")
    init_density_concentration = sp.Symbol("init_density_concentration")
    init_velocity_fluid = sp.symbols("init_velocity_fluid_:3")
    #init_velocity_concentration = sp.symbols("init_velocity_concentration_:3")
    pdfs_inter_fluid = sp.symbols("pdfs_inter_fluid:" + str(stencil_fluid.Q))
    pdfs_inter_concentration = sp.symbols("pdfs_inter_concentration:" + str(stencil_concentration.Q))
    rho_0 = sp.Symbol("rho_0")
    T0 = sp.Symbol("T0")
    alpha = sp.Symbol("alpha")
    gravity_LBM = sp.Symbol("gravityLB")
    Sv = sp.Symbol("Sv")
    Sq = sp.Symbol("Sq")
    omega_f = sp.Symbol("omega_f")
    omega_c = sp.Symbol("omega_c")



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


    # Solid collision variant
    SC = int(config_tokens[1])
    qk = sp.Symbol("qk")
    qe = sp.Symbol("qe")

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

    pdfs_energy, pdfs_energy_tmp, energy_field = ps.fields(
        f"pdfs_energy({stencil_energy.Q}), pdfs_energy_tmp({stencil_energy.Q}), energy_field({1}): {data_type}[3D]",
        layout=layout,
    )

    # particle related fields (considering for all the particle, i.e: MaxParticlesPerCell

    particle_velocities, particle_forces, Bs = ps.fields(
        f"particle_v({MaxParticlesPerCell * stencil_fluid.D}), particle_f({MaxParticlesPerCell * stencil_fluid.D}), Bs({MaxParticlesPerCell}): {data_type}[3D]",
        layout=layout,
    )
    particle_temperatures = ps.fields(f"particle_t({MaxParticlesPerCell}) :{data_type}[3D]", layout=layout)
    particle_energies = ps.fields(f"particle_energy({MaxParticlesPerCell}) :{data_type}[3D]", layout=layout)

    # Solid fraction field
    B = ps.fields(f"b({1}): {data_type}[3D]", layout=layout)


    #force_concentration_on_fluid = sp.Matrix([0, (rho_0)*alpha*(concentration_field.center - T0)*gravity_LBM,0])
    force_concentration_on_fluid = sp.Matrix([0, 0,(rho_0)*alpha*(concentration_field.center - T0)*gravity_LBM])

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

    # Energy LBM optimisation
    lbm_energy_opt = LBMOptimisation(
        cse_global=False,
        symbolic_field=pdfs_energy,
        symbolic_temporary_field=pdfs_energy_tmp,
        field_layout=layout,
    )

    # Fluid LBM config
    lbm_fluid_config = LBMConfig(
        stencil=stencil_fluid,
        method=Method.SRT,
        relaxation_rate=omega_f,
        output={"velocity": velocity_field},
        force= force_concentration_on_fluid,
        force_model=ForceModel.GUO,
        compressible=True,
    )

    # Fluid PSM config
    psm_config_F = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        SC=SC,
        MaxParticlesPerCell=MaxParticlesPerCell,
        individual_fraction_field=Bs,
        particle_force_field=particle_forces,
    )

    psm_fluid_config = LBMConfig(
        stencil=stencil_fluid,
        method=Method.SRT,
        relaxation_rate=omega_f,
        output={"velocity": velocity_field},
        force= force_concentration_on_fluid,
        force_model=ForceModel.LUO,
        compressible=False,
        psm_config=psm_config_F,
    )

    # Concentration LBM config
    lbm_concentration_config = LBMConfig(
        stencil=stencil_concentration,
        method=Method.MRT,
        relaxation_rates=[1,qk,qk,qk,qe,qe,qe],  # MRT 2
        #relaxation_rate=omega_c,
        velocity_input=velocity_field,
        output={"density": concentration_field},
        compressible=True,
        zero_centered=False
    )

    # Concentration PSM config
    psm_config_C = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        SC=int(4),
        MaxParticlesPerCell=MaxParticlesPerCell,
        individual_fraction_field=Bs,
        particle_force_field=None,
        particle_temperature_field=particle_temperatures,
    )


    psm_concentration_config = LBMConfig(
        stencil=stencil_concentration,
        method=Method.SRT,
        relaxation_rate=omega_c,
        velocity_input=velocity_field,
        output={"density": concentration_field},
        compressible=True,
        psm_config=psm_config_C,
        continuous_equilibrium=False,
        zero_centered=False,
    )

    ## for CHT
    omega_p = sp.Symbol("omega_p")  # for the particles in CHT
    rho_f = sp.Symbol("rho_f")
    rho_s = sp.Symbol("rho_s")
    omegaT_f = sp.Symbol("omegaT_f")
    omegaT_s = sp.Symbol("omegaT_s")
    Cp_f = sp.Symbol("Cp_f")
    Cp_s = sp.Symbol("Cp_s")

    # Energy PSM config
    psm_config_E = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        SC=int(5),
        MaxParticlesPerCell=MaxParticlesPerCell,
        individual_fraction_field=Bs,
        particle_force_field=None,
        particle_temperature_field=particle_temperatures,
        particle_density=rho_s,
        particle_specific_heat=Cp_s,
        solid_relaxation_rate=omegaT_s,
        energy_field=energy_field,
        temperature_field=concentration_field,
    )


    psm_energy_config = LBMConfig(
        stencil=stencil_energy,
        method=Method.SRT,
        relaxation_rate=omegaT_f,  # omega_f will be used for the fluid and omega_p will be used for the solid particles
        velocity_input=velocity_field,
        output={"density": energy_field},
        compressible=True,
        psm_config=psm_config_E,
        continuous_equilibrium=False,
        zero_centered=False,
    )


    if config_tokens[0] == "srt-smagorinsky" or config_tokens[0] == "trt-smagorinsky":
        lbm_fluid_config.smagorinsky = True

    # =====================
    # Generate method
    # =====================

    method_fluid = create_lb_method(lbm_config=psm_fluid_config)
    method_concentration = create_lb_method(lbm_config=psm_concentration_config)
    method_energy = create_lb_method(lbm_config=psm_energy_config)
    init_velocity = sp.symbols("init_velocity_:3")
    pdfs_fluid_setter = macroscopic_values_setter(
        method_fluid, density=init_density_fluid, velocity=velocity_field.center_vector, pdfs=pdfs_fluid.center_vector
    )

    pdfs_concentration_setter = macroscopic_values_setter(
        method_concentration, density=concentration_field.center, velocity= velocity_field.center_vector,pdfs=pdfs_concentration.center_vector
    )


    pdfs_energy_setter = macroscopic_values_setter(
        method_energy, density=energy_field.center, velocity= velocity_field.center_vector,pdfs=pdfs_energy.center_vector
    )

    # Use average velocity of all intersecting particles when setting PDFs (mandatory for SC=3)
    rhs = []
    for i, sub_exp in enumerate(pdfs_fluid_setter.subexpressions[-3:]):
        if(len(sub_exp.rhs.args) > 0):
            for summand in (sub_exp.rhs.args):
                rhs.append(summand * (1.0 - B.center))
        else:
            rhs.append(sub_exp.rhs * (1.0 - B.center))
        for p in range(2):
            rhs.append(particle_velocities(p * stencil_fluid.D + i) * Bs.center(p))
        pdfs_fluid_setter.subexpressions.remove(sub_exp)
        pdfs_fluid_setter.subexpressions.append(Assignment(sub_exp.lhs, Add(*rhs)))
        rhs = []


    # Use average temperature of all intersecting particles when setting PDFs (mandatory for SC=3)

    sub_exp_con = pdfs_concentration_setter.subexpressions[0]
    rhs_con = []
    rhs_con.append(sub_exp_con.rhs * (1.0 - B.center))
    for p in range(2):
        rhs_con.append(particle_temperatures(p) * Bs.center(p))
    pdfs_concentration_setter.subexpressions.remove(sub_exp_con)
    pdfs_concentration_setter.subexpressions.append(Assignment(sub_exp_con.lhs, Add(*rhs_con)))
    print("con setter after manip ", pdfs_concentration_setter.subexpressions)


    ## for energy
    sub_exp_energy = pdfs_energy_setter.subexpressions[0]
    rhs_energy = []
    rhs_energy.append((rho_f*Cp_f*omegaT_f*(1-B.center)*concentration_field.center)/((1 - B.center)*omegaT_f + B.center*omegaT_s))
    for p in range(2):
        rhs_energy.append((particle_temperatures(p) * Bs.center(p) * rho_s * Cp_s * omegaT_s)/((1 - B.center)*omegaT_f + B.center*omegaT_s))
    pdfs_energy_setter.subexpressions.remove(sub_exp_energy)
    pdfs_energy_setter.subexpressions.append(Assignment(sub_exp_energy.lhs, Add(*rhs_energy)))

    #print("energy setter after manip ", pdfs_energy_setter.subexpressions)

    # specify the target


    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU

    node_collection_fluid = create_psm_update_rule(lbm_config=psm_fluid_config, lbm_optimisation=lbm_fluid_opt)
    node_collection_concentration = create_psm_update_rule(lbm_config=psm_concentration_config, lbm_optimisation=lbm_concentration_opt)
    node_collection_energy = create_psm_update_rule(lbm_config=psm_energy_config, lbm_optimisation=lbm_energy_opt)

    ## defining custom pystencils kernel that computes temperature from rho_cp_T

    rho_cp_eff = sp.Symbol("rho_cp_eff")
    @ps.kernel
    def compute_temperature_field():

        rho_cp_eff = ((1.0 - B.center)* method_fluid.conserved_quantity_computation.density_symbol *Cp_f*omegaT_f + B.center*rho_s*Cp_s*omegaT_s)/((1-B.center)*omegaT_f + B.center*omegaT_s)
        concentration_field.center @= energy_field.center/rho_cp_eff
    compute_temperature_field_ac = ps.AssignmentCollection(
            compute_temperature_field
    )
    generate_sweep(ctx, "compute_temperature_field", compute_temperature_field_ac,target=target)


    assignments = []
    assignments.append(method_energy.conserved_quantity_computation.equilibrium_input_equations_from_pdfs(pdfs_energy.center_vector))
    for k in range(MaxParticlesPerCell):
        rho_cp_eff = ((1.0 - B.center)* method_fluid.conserved_quantity_computation.density_symbol *Cp_f*omegaT_f + B.center*rho_s*Cp_s*omegaT_s)/((1-B.center)*omegaT_f + B.center*omegaT_s)
        condition = sp.Piecewise(
            (method_energy.conserved_quantity_computation.density_symbol/rho_cp_eff , B.center > 0),
            (0, True)
        )
        assignments.append(ps.Assignment(particle_temperatures.center(k), condition))

    ac = ps.AssignmentCollection(assignments)
    generate_sweep(ctx, "compute_temperature_field_particle", ac)


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
        "PSMFluidSweep",
        node_collection_fluid,
        field_swaps=[(pdfs_fluid, pdfs_fluid_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "PSMConcentrationSweep",
        node_collection_concentration,
        field_swaps=[(pdfs_concentration, pdfs_concentration_tmp)],
        target=target,
    )

    generate_sweep(
        ctx,
        "PSMEnergySweep",
        node_collection_energy,
        field_swaps=[(pdfs_energy, pdfs_energy_tmp)],
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
    generate_sweep(
        ctx,
        "PSMFluidSweepSplit",
        node_collection_fluid,
        field_swaps=[(pdfs_fluid, pdfs_fluid_tmp)],
        target=target,
        inner_outer_split=True,
    )

    generate_sweep(
        ctx,
        "PSMConcentrationSweepSplit",
        node_collection_concentration,
        field_swaps=[(pdfs_concentration, pdfs_concentration_tmp)],
        target=target,
        inner_outer_split=True,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PackInfoFluid",
        create_lb_update_rule(lbm_config=psm_fluid_config, lbm_optimisation=lbm_fluid_opt),
        target=target,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PackInfoConcentration",
        create_lb_update_rule(lbm_config=psm_concentration_config, lbm_optimisation=lbm_concentration_opt),
        target=target,
    )

    generate_pack_info_from_kernel(
        ctx,
        "PackInfoEnergy",
        create_lb_update_rule(lbm_config=psm_energy_config, lbm_optimisation=lbm_energy_opt),
        target=target,
    )

    generate_sweep(ctx, "InitializeFluidDomain", pdfs_fluid_setter, target=target)
    generate_sweep(ctx, "InitializeConcentrationDomain", pdfs_concentration_setter, target=target)
    generate_sweep(ctx, "InitializeEnergyDomain", pdfs_energy_setter, target=target)

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

    generate_boundary(
        ctx,
        "BC_Fluid_Outflow",
        ExtrapolationOutflow((0,0,1),method_fluid),
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
        target=target, )



    # energy boundary conditions

    dirichlet_bc_dynamic = DiffusionDirichlet(lambda *args: None, velocity_field, data_type=data_type)
    diffusion_data_handler = DiffusionDirichletAdditionalDataHandler(stencil_energy, dirichlet_bc_dynamic)
    generate_boundary(ctx, 'BC_energy_DiffusionDirichlet_dynamic', dirichlet_bc_dynamic, method_energy,
                      additional_data_handler=diffusion_data_handler,
                      target=target, streaming_pattern='pull', data_type=data_type)

    bc_density_energy = sp.Symbol("bc_density_energy")
    generate_boundary(
        ctx,
        "BC_energy_DiffusionDirichlet_static",
        DiffusionDirichlet(bc_density_energy),
        method_energy,
        field_name=pdfs_energy.name,
        streaming_pattern="pull",
        target=target,
    )


    generate_boundary(
        ctx,
        "BC_Energy_Neumann",
        NeumannByCopy(stencil_energy),
        method_energy,
        field_name=pdfs_energy.name,
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


    stencil_typedefs = {"Stencil_Fluid_T": stencil_fluid, "CommunicationStencil_Fluid_T": stencil_fluid, "Stencil_Concentration_T": stencil_concentration, "CommunicationStencil_Concentration_T": stencil_concentration
                        , "Stencil_Energy_T":stencil_energy, "CommunicationStencil_Energy_T":stencil_energy}
    field_typedefs = {
        "PdfField_fluid_T": pdfs_fluid,
        "DensityField_fluid_T": density_field,
        "VelocityField_fluid_T": velocity_field,
        "PdfField_concentration_T": pdfs_concentration,
        "DensityField_concentration_T": concentration_field,
        "PdfField_energy_T": pdfs_energy,
        "DensityField_energy_T": energy_field,
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

    pdfs_energy_getter = macroscopic_values_getter(
       method_energy, density=energy_field, velocity=None,pdfs=pdfs_energy.center_vector
    )

    generate_sweep(ctx, "FluidMacroSetter", pdfs_fluid_setter)
    generate_sweep(ctx, "FluidMacroGetter", pdfs_fluid_getter)
    generate_sweep(ctx, "ConcentrationMacroSetter", pdfs_concentration_setter)
    generate_sweep(ctx, "ConcentrationMacroGetter", pdfs_concentration_getter)
    generate_sweep(ctx, "EnergyMacroSetter", pdfs_energy_setter)
    generate_sweep(ctx, "EnergyMacroGetter", pdfs_energy_getter)