from pystencils.field import fields
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Stencil, Method
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow
from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler, OutflowAdditionalDataHandler
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header
from lbmpy_walberla import generate_boundary, generate_lb_pack_info

import sympy as sp

stencil = LBStencil(Stencil.D2Q9)

pdfs, pdfs_tmp = fields(f"pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): double[{stencil.D}D]", layout='fzyx')
velocity_field, density_field = fields(f"velocity({stencil.D}), density(1) : double[{stencil.D}D]", layout='fzyx')
omega = sp.Symbol("omega")
u_max = sp.Symbol("u_max")

output = {
    'density': density_field,
    'velocity': velocity_field
}

lbm_config = LBMConfig(method=Method.CUMULANT, stencil=stencil, relaxation_rate=omega,
                       galilean_correction=stencil.Q == 27, field_name='pdfs', output=output)

lbm_opt = LBMOptimisation(symbolic_field=pdfs, symbolic_temporary_field=pdfs_tmp,
                          cse_global=False, cse_pdfs=False)

method = create_lb_method(lbm_config=lbm_config)

# getter & setter
setter_assignments = macroscopic_values_setter(method, velocity=velocity_field.center_vector,
                                               pdfs=pdfs, density=1.0)

update_rule = create_lb_update_rule(lb_method=method, lbm_config=lbm_config, lbm_optimisation=lbm_opt)

stencil_typedefs = {'Stencil_T': stencil}
field_typedefs = {'PdfField_T': pdfs,
                  'VelocityField_T': velocity_field,
                  'ScalarField_T': density_field}

with CodeGeneration() as ctx:
    # sweeps
    generate_sweep(ctx, 'GeneratedOutflowBC_Sweep', update_rule, field_swaps=[(pdfs, pdfs_tmp)])
    generate_sweep(ctx, 'GeneratedOutflowBC_MacroSetter', setter_assignments)

    # boundaries
    ubb_dynamic = UBB(lambda *args: None, dim=stencil.D)
    ubb_data_handler = UBBAdditionalDataHandler(stencil, ubb_dynamic)

    if stencil.D == 2:
        ubb_static = UBB([sp.Symbol("u_max"), 0])
    else:
        ubb_static = UBB([sp.Symbol("u_max"), 0, 0])

    outflow = ExtrapolationOutflow(stencil[4], method)
    outflow_data_handler = OutflowAdditionalDataHandler(stencil, outflow)

    # Dynamic UBB which is used to produce a specific velocity profile at the inflow.
    # Note that the additional data handler is needed for that kind of boundary.
    generate_boundary(ctx, 'GeneratedOutflowBC_Dynamic_UBB', ubb_dynamic, method,
                      additional_data_handler=ubb_data_handler)

    # Static UBB which is used to apply a certain velocity u_max at the upper wall in x-direction
    generate_boundary(ctx, 'GeneratedOutflowBC_Static_UBB', ubb_static, method)

    generate_boundary(ctx, 'GeneratedOutflowBC_NoSlip', NoSlip(), method)

    generate_boundary(ctx, 'GeneratedOutflowBC_Outflow', outflow, method,
                      additional_data_handler=outflow_data_handler)

    # communication
    generate_lb_pack_info(ctx, 'GeneratedOutflowBC_PackInfo', stencil, pdfs)

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, "GeneratedOutflowBC.h",
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
