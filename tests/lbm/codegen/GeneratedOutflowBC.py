from pystencils.field import fields
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.stencils import get_stencil
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule
from lbmpy.boundaries import NoSlip, UBB, ExtrapolationOutflow
from lbmpy_walberla.additional_data_handler import UBBAdditionalDataHandler, OutflowAdditionalDataHandler
from pystencils_walberla import CodeGeneration, generate_sweep
from lbmpy_walberla import RefinementScaling, generate_boundary, generate_lb_pack_info

import sympy as sp

stencil = get_stencil("D2Q9")
q = len(stencil)
dim = len(stencil[0])

pdfs, pdfs_tmp = fields(f"pdfs({q}), pdfs_tmp({q}): double[{dim}D]", layout='fzyx')
velocity_field, density_field = fields(f"velocity({dim}), density(1) : double[{dim}D]", layout='fzyx')
omega = sp.Symbol("omega")
u_max = sp.Symbol("u_max")

output = {
    'density': density_field,
    'velocity': velocity_field
}

options = {'method': 'cumulant',
           'stencil': stencil,
           'relaxation_rate': omega,
           'galilean_correction': len(stencil) == 27,
           'field_name': 'pdfs',
           'output': output,
           'optimization': {'symbolic_field': pdfs,
                            'symbolic_temporary_field': pdfs_tmp,
                            'cse_global': False,
                            'cse_pdfs': False}}

method = create_lb_method(**options)

# getter & setter
setter_assignments = macroscopic_values_setter(method, velocity=velocity_field.center_vector,
                                               pdfs=pdfs, density=1.0)

# opt = {'instruction_set': 'sse', 'assume_aligned': True, 'nontemporal': False, 'assume_inner_stride_one': True}

update_rule = create_lb_update_rule(lb_method=method, **options)

info_header = f"""
using namespace walberla;
#include "stencil/D{dim}Q{q}.h"
using Stencil_T = walberla::stencil::D{dim}Q{q};
using PdfField_T = GhostLayerField<real_t, {q}>;
using VelocityField_T = GhostLayerField<real_t, {dim}>;
using ScalarField_T = GhostLayerField<real_t, 1>;
    """

stencil = method.stencil

with CodeGeneration() as ctx:
    # sweeps
    generate_sweep(ctx, 'GeneratedOutflowBC_Sweep', update_rule, field_swaps=[(pdfs, pdfs_tmp)])
    generate_sweep(ctx, 'GeneratedOutflowBC_MacroSetter', setter_assignments)

    # boundaries
    ubb_dynamic = UBB(lambda *args: None, dim=dim)
    ubb_data_handler = UBBAdditionalDataHandler(stencil, ubb_dynamic)

    if dim == 2:
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
    ctx.write_file("GeneratedOutflowBC_InfoHeader.h", info_header)
