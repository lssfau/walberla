import sympy as sp
from pystencils import Target

from pystencils.field import fields
from pystencils.simp import insert_aliases, insert_constants

from lbmpy import LBStencil, LBMConfig, LBMOptimisation
from lbmpy.boundaries.boundaryconditions import (ExtrapolationOutflow, UBB, QuadraticBounceBack,
                                                 FreeSlip, NoSlip, FixedDensity)
from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.enums import Method, Stencil, SubgridScaleModel

from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

info_header = """
#pragma once
#include <map>
std::map<std::string, std::string> infoMap{{{{"stencil", "{stencil}"}},
                                           {{"streamingPattern", "{streaming_pattern}"}}, 
                                           {{"collisionOperator", "{collision_operator}"}}}};
"""

omega = sp.symbols("omega")
inlet_velocity = sp.symbols("u_x")
max_threads = 256

with CodeGeneration() as ctx:
    target = Target.GPU if ctx.gpu else Target.CPU
    sweep_params = {'block_size': (128, 1, 1)} if ctx.gpu else {}

    # The application is meant to be compiled with double precision. For single precision, the pdf_dtype can be switched
    # to float32. In this way the calculations are still performed in double precision while the application can profit
    # from enhanced performance due to the lower precision of the PDF field
    dtype = 'float64' if ctx.double_accuracy else 'float32'
    pdf_dtype = 'float64'

    stencil = LBStencil(Stencil.D3Q27)
    q = stencil.Q
    dim = stencil.D

    streaming_pattern = 'esopull'

    pdfs = fields(f"pdfs_src({stencil.Q}): {pdf_dtype}[3D]", layout='fzyx')
    velocity_field, density_field = fields(f"velocity({dim}), density(1) : {dtype}[{dim}D]", layout='fzyx')
    omega_field = fields(f"rr(1) : {dtype}[{dim}D]", layout='fzyx')

    macroscopic_fields = {'density': density_field, 'velocity': velocity_field}

    method_enum = Method.CUMULANT
    lbm_config = LBMConfig(
        method=method_enum,
        stencil=stencil,
        relaxation_rate=omega_field.center,
        compressible=True,
        # subgrid_scale_model=SubgridScaleModel.QR,
        fourth_order_correction=0.01 if method_enum == Method.CUMULANT and stencil.Q == 27 else False,
        field_name='pdfs',
        streaming_pattern=streaming_pattern,
    )

    lbm_opt = LBMOptimisation(cse_global=False, cse_pdfs=False, field_layout="fzyx",
                              symbolic_field=pdfs)

    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    if lbm_config.method == Method.CUMULANT:
        collision_rule = insert_constants(collision_rule)
        collision_rule = insert_aliases(collision_rule)
    lb_method = collision_rule.method

    freeslip = lbm_boundary_generator("FreeSlip", flag_uid="FreeSlip", boundary_object=FreeSlip(stencil),
                                      field_data_type=pdf_dtype)
    no_slip = lbm_boundary_generator(class_name='NoSlip', flag_uid='NoSlip',
                                     boundary_object=NoSlip(), field_data_type=pdf_dtype)

    quadratic_bounce_back = QuadraticBounceBack(omega, calculate_force_on_boundary=True, data_type=dtype)
    no_slip_interpolated = lbm_boundary_generator(class_name='Obstacle', flag_uid='Obstacle',
                                                  boundary_object=quadratic_bounce_back,
                                                  field_data_type=pdf_dtype)

    ubb = lbm_boundary_generator(class_name='UBB', flag_uid='UBB',
                                 boundary_object=UBB((inlet_velocity, 0.0, 0.0), density=1.0, data_type=dtype),
                                 field_data_type=pdf_dtype)

    outflow_boundary = ExtrapolationOutflow(stencil[4], lb_method, data_type=pdf_dtype)
    outflow = lbm_boundary_generator(class_name='Outflow', flag_uid='Outflow',
                                     boundary_object=outflow_boundary,
                                     field_data_type=pdf_dtype)

    fixed_density = lbm_boundary_generator(class_name='FixedDensity', flag_uid='FixedDensity',
                                           boundary_object=FixedDensity(1.0),
                                           field_data_type=pdf_dtype)

    generate_lbm_package(ctx, name="FlowAroundSphere", collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=True, boundaries=[freeslip, no_slip, no_slip_interpolated,
                                                      ubb, outflow, fixed_density],
                         macroscopic_fields=macroscopic_fields, gpu_indexing_params=sweep_params,
                         target=target, data_type=dtype, pdfs_data_type=pdf_dtype,
                         max_threads=max_threads)

    field_typedefs = {'VelocityField_T': velocity_field,
                      'ScalarField_T': density_field}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'FlowAroundSphereInfoHeader',
                         field_typedefs=field_typedefs)

    infoHeaderParams = {
        'stencil': stencil.name.lower(),
        'streaming_pattern': streaming_pattern,
        'collision_operator': lbm_config.method.name.lower(),
    }

    ctx.write_file("FlowAroundSphereStaticDefines.h", info_header.format(**infoHeaderParams))
