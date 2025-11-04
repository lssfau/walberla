import sympy as sp
import pystencils as ps

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.partially_saturated_cells import PSMConfig

from lbmpy.boundaries import NoSlip, UBB, FixedDensity, FreeSlip
from lbmpy.creationfunctions import (
    create_lb_collision_rule,
    create_lb_method,
    create_psm_update_rule,
)

from lbmpy.macroscopic_value_kernels import (
    macroscopic_values_setter,
    macroscopic_values_getter,
)

from pystencils_walberla import (
    CodeGeneration,
    generate_info_header,
    generate_sweep,
)

from lbmpy_walberla import (
    generate_lbm_package,
    lbm_boundary_generator,
)

info_header = """
const char * infoStencil = "{stencil}";
const char * infoStreamingPattern = "{streaming_pattern}";
const char * infoCollisionSetup = "{collision_setup}";
const bool infoCseGlobal = {cse_global};
const bool infoCsePdfs = {cse_pdfs};
"""

with CodeGeneration() as ctx:
    # =====================
    # Parameters
    # =====================
    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)
    assert stencil.D == 3
    omega = sp.Symbol("omega")
    init_density = sp.Symbol("init_density")
    init_velocity = sp.symbols("init_velocity_:3")
    pdfs_inter = sp.symbols("pdfs_inter:" + str(stencil.Q))
    layout = "fzyx"
    config_tokens = ctx.config.split("_")
    methods = {
        "srt": Method.SRT,
        "trt": Method.TRT,
        "mrt": Method.MRT,
        "cumulant": Method.CUMULANT,
        "srt-smagorinsky": Method.SRT,
        "trt-smagorinsky": Method.TRT,
    }
    # Solid collision variant
    SC = int(config_tokens[1][2])
    # Force model
    kwargs = {}
    if len(config_tokens) == 3:
        MaxParticlesPerCell = int(config_tokens[2])
    elif len(config_tokens) == 4 and config_tokens[2] == "force":
        kwargs["force"] = sp.symbols("F_:3")
        kwargs["force_model"] = ForceModel.LUO
        MaxParticlesPerCell = int(config_tokens[3])
    else:
        raise ValueError("Wrong configuration!")

    # =====================
    # Fields
    # =====================
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

    # =====================
    # Lbmpy settings
    # =====================
    psm_config = PSMConfig(
        fraction_field=B,
        object_velocity_field=particle_velocities,
        solid_collision=SC,
        max_particles_per_cell=MaxParticlesPerCell,
        individual_fraction_field=Bs,
        object_force_field=particle_forces,
    )

    lbm_config = LBMConfig(
        stencil=stencil,
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        compressible=True,
        psm_config=psm_config,
        **kwargs,
    )

    lbm_opt = LBMOptimisation(
        symbolic_field=pdfs,
        symbolic_temporary_field=pdfs_tmp,
        field_layout=layout,
        cse_global=True,
    )
    if config_tokens[0] == "srt-smagorinsky" or config_tokens[0] == "trt-smagorinsky":
        lbm_config.smagorinsky = 0.1

    # =====================
    # Generate PSM update sweep
    # =====================
    method = create_lb_method(lbm_config=lbm_config)
    node_collection = create_psm_update_rule(lbm_config, lbm_opt)
    setter_assignments = macroscopic_values_setter(
        method, init_density, init_velocity, pdfs.center_vector, psm_config=psm_config
    )

    getter_assignments = macroscopic_values_getter(
        method,
        density=density_field,
        velocity=velocity_field.center_vector,
        pdfs=pdfs.center_vector,
        psm_config=psm_config,
    )

    # Getter & setter to compute moments from pdfs
    generate_sweep(ctx, "PSM_MacroGetter", getter_assignments, target=target)
    generate_sweep(ctx, "PSM_MacroSetter", setter_assignments, target=target)

    generate_sweep(
        ctx,
        "PSMSweep",
        node_collection,
        field_swaps=[(pdfs, pdfs_tmp)],
        target=target,
        inner_outer_split=True,
    )

    # =====================
    # Generate LBM package
    # =====================
    # Reinitialize the LBM config without PSM
    lbm_config = LBMConfig(
        stencil=stencil,
        method=methods[config_tokens[0]],
        relaxation_rate=omega,
        compressible=True,
        **kwargs,
    )
    if config_tokens[0] == "srt-smagorinsky" or config_tokens[0] == "trt-smagorinsky":
        lbm_config.smagorinsky = 0.1

    # Used for setting PDFs according to the velocity field
    method = create_lb_method(lbm_config=lbm_config)
    setter_assignments_vel_field = macroscopic_values_setter(
        method,
        density=1.0,
        velocity=velocity_field.center_vector,
        pdfs=pdfs.center_vector,
    )

    no_slip = lbm_boundary_generator(
        class_name="LBM_NoSlip", flag_uid="NoSlip", boundary_object=NoSlip()
    )
    free_slip = lbm_boundary_generator(
        class_name="LBM_FreeSlip",
        flag_uid="FreeSlip",
        boundary_object=FreeSlip(stencil=stencil),
    )
    bc_velocity = sp.symbols("bc_velocity_:3")
    ubb = lbm_boundary_generator(
        class_name="LBM_UBB",
        flag_uid="UBB",
        boundary_object=UBB(bc_velocity, data_type=data_type),
    )
    bc_density = sp.Symbol("bc_density")
    bc_density_2 = sp.Symbol("bc_density_2")
    fixed_density = lbm_boundary_generator(
        class_name="LBM_FixedDensity",
        flag_uid="FixedDensity",
        boundary_object=FixedDensity(bc_density),
    )
    fixed_density_2 = lbm_boundary_generator(
        class_name="LBM_FixedDensity_2",
        flag_uid="FixedDensity2",
        boundary_object=FixedDensity(bc_density_2),
    )

    macroscopic_fields = {"density": density_field, "velocity": velocity_field}

    generate_lbm_package(
        ctx,
        name="LBM",
        collision_rule=create_lb_collision_rule(
            lbm_config=lbm_config, lbm_optimisation=lbm_opt
        ),
        lbm_config=lbm_config,
        lbm_optimisation=lbm_opt,
        boundaries=[
            no_slip,
            free_slip,
            ubb,
            fixed_density,
            fixed_density_2,
        ],
        macroscopic_fields=macroscopic_fields,
        target=target,
    )
    generate_sweep(ctx, "LBM_MacroSetter", setter_assignments_vel_field, target=target)

    stencil_typedefs = {"Stencil_T": stencil, "CommunicationStencil_T": stencil}
    field_typedefs = {
        "PdfField_T": pdfs,
        "DensityField_T": density_field,
        "VelocityField_T": velocity_field,
    }

    infoHeaderParams = {
        "stencil": stencil.name,
        "streaming_pattern": lbm_config.streaming_pattern,
        "collision_setup": config_tokens[0],
        "cse_global": int(lbm_opt.cse_global),
        "cse_pdfs": int(lbm_opt.cse_pdfs),
    }

    generate_info_header(
        ctx,
        "LBM_InfoHeader",
        stencil_typedefs=stencil_typedefs,
        field_typedefs=field_typedefs,
        additional_code=info_header.format(**infoHeaderParams),
    )
