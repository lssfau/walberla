import sympy as sp
import pystencils as ps

from lbmpy.enums import SubgridScaleModel
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil, ForceModel
from lbmpy.flow_statistics import welford_assignments
from lbmpy.relaxationrates import lattice_viscosity_from_relaxation_rate

from lbmpy.creationfunctions import create_lb_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

from lbmpy.boundaries import NoSlip, FreeSlip, WallFunctionBounce
from lbmpy.boundaries.wall_function_models import SpaldingsLaw, MoninObukhovSimilarityTheory
from lbmpy.utils import frobenius_norm, second_order_moment_tensor

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel
from lbmpy_walberla import generate_boundary
from lbmpy_walberla.additional_data_handler import WFBAdditionalDataHandler

#   =====================
#      Code Generation
#   =====================

info_header = """
#ifndef TURBULENTCHANNEL_INCLUDES
#define TURBULENTCHANNEL_INCLUDES

#include <stencil/D{d}Q{q}.h>

#include "TurbulentChannel_Sweep.h"
#include "TurbulentChannel_PackInfo.h"
#include "TurbulentChannel_Setter.h"
#include "TurbulentChannel_Welford.h"
#include "TurbulentChannel_Welford_TKE_SGS.h"
#include "TurbulentChannel_TKE_SGS_Writer.h"

#include "TurbulentChannel_NoSlip.h"
#include "TurbulentChannel_FreeSlip_top.h"
#include "TurbulentChannel_WFB_top.h"
#include "TurbulentChannel_WFB_bottom.h"

{config_definitions}

namespace walberla {{
    namespace codegen {{
        using Stencil_T = walberla::stencil::D{d}Q{q};
        static constexpr uint_t flowAxis = {flow_axis};
        static constexpr uint_t wallAxis = {wall_axis};
        
        static constexpr field::Layout layout = field::{layout};
    }}
}}

#endif // TURBULENTCHANNEL_INCLUDES
"""


def check_axis(flow_axis, wall_axis):
    assert flow_axis != wall_axis
    assert flow_axis < 3
    assert wall_axis < 3


with CodeGeneration() as ctx:
    ctx.optimize_for_localhost = False

    flow_axis = 0
    wall_axis = 1

    check_axis(flow_axis=flow_axis, wall_axis=wall_axis)

    #   ========================
    #      General Parameters
    #   ========================

    config_tokens = ctx.config.split('_')

    assert len(config_tokens) >= 4
    stencil_str = config_tokens[0]
    wfb_weights_str = config_tokens[1]
    velocity_str = config_tokens[2]
    sgs_str = config_tokens[3]

    target = ps.Target.CPU

    data_type = "float64" if ctx.double_accuracy else "float32"
    if stencil_str == 'd3q19':
        stencil = LBStencil(Stencil.D3Q19)
    elif stencil_str == 'd3q27':
        stencil = LBStencil(Stencil.D3Q27)
    else:
        raise ValueError("Only D3Q27 and D3Q19 stencil are supported at the moment")
    omega = sp.Symbol('omega')

    F_x = sp.Symbol('F_x')
    force_vector = [0] * 3
    force_vector[flow_axis] = F_x

    layout = 'fzyx'

    normal_direction_top = [0] * 3
    normal_direction_top[wall_axis] = -1
    normal_direction_top = tuple(normal_direction_top)

    normal_direction_bottom = [0] * 3
    normal_direction_bottom[wall_axis] = 1
    normal_direction_bottom = tuple(normal_direction_bottom)

    #   PDF Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]', layout=layout)

    #   Output Fields
    omega_field = ps.fields(f"omega_out: {data_type}[{stencil.D}D]", layout=layout)
    sgs_tke = ps.fields(f"sgs_tke: {data_type}[{stencil.D}D]", layout=layout)
    mean_sgs_tke = ps.fields(f"mean_sgs_tke: {data_type}[{stencil.D}D]", layout=layout)
    velocity = ps.fields(f"velocity({stencil.D}): {data_type}[{stencil.D}D]", layout=layout)
    mean_velocity = ps.fields(f"mean_velocity({stencil.D}): {data_type}[{stencil.D}D]", layout=layout)
    sum_of_squares_field = ps.fields(f"sum_of_squares_field({stencil.D**2}): {data_type}[{stencil.D}D]", layout=layout)

    # LBM Optimisation
    lbm_opt = LBMOptimisation(cse_global=True,
                              symbolic_field=pdfs,
                              symbolic_temporary_field=pdfs_tmp,
                              field_layout=layout)

    #   ==================
    #      Method Setup
    #   ==================

    lbm_config = LBMConfig(stencil=stencil,
                           method=Method.CUMULANT,
                           force_model=ForceModel.GUO,
                           force=tuple(force_vector),
                           relaxation_rate=omega,
                           compressible=True,
                           output={'velocity': velocity})

    if stencil_str == 'd3q27':
        lbm_config.galilean_correction = True
        lbm_config.fourth_order_correction = True

    if sgs_str == 'smago':
        lbm_config.subgrid_scale_model = SubgridScaleModel.SMAGORINSKY
    elif sgs_str == 'qr':
        lbm_config.subgrid_scale_model = SubgridScaleModel.QR

    if not sgs_str == 'none':
        lbm_config.omega_output_field = omega_field

    update_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    lbm_method = update_rule.method

    #   ========================
    #      PDF Initialization
    #   ========================

    initial_rho = sp.Symbol('rho_0')

    pdfs_setter = macroscopic_values_setter(lbm_method,
                                            initial_rho,
                                            velocity.center_vector,
                                            pdfs.center_vector)

    #   LBM Sweep
    generate_sweep(ctx, "TurbulentChannel_Sweep", update_rule, field_swaps=[(pdfs, pdfs_tmp)], target=target)

    #   Pack Info
    generate_pack_info_from_kernel(ctx, "TurbulentChannel_PackInfo", update_rule, target=target)

    #   Macroscopic Values Setter
    generate_sweep(ctx, "TurbulentChannel_Setter", pdfs_setter, target=target, ghost_layers_to_include=1)

    #   Welford update
    welford_update = welford_assignments(field=velocity, mean_field=mean_velocity,
                                         sum_of_squares_field=sum_of_squares_field)
    generate_sweep(ctx, "TurbulentChannel_Welford", welford_update, target=target)

    tke_welford_update = welford_assignments(field=sgs_tke, mean_field=mean_sgs_tke)
    generate_sweep(ctx, "TurbulentChannel_Welford_TKE_SGS", tke_welford_update, target=target)

    # subgrid TKE output
    @ps.kernel
    def tke_sgs_writer():
        f_neq = sp.Matrix(pdfs.center_vector) - lbm_method.get_equilibrium_terms()
        rho = lbm_method.conserved_quantity_computation.density_symbol
        strain_rate = frobenius_norm(-3 * omega_field.center / (2 * rho) * second_order_moment_tensor(f_neq, lbm_method.stencil))
        eddy_viscosity = lattice_viscosity_from_relaxation_rate(omega_field.center) - lattice_viscosity_from_relaxation_rate(omega)

        sgs_tke.center @= (eddy_viscosity * strain_rate**2)**(2.0/3.0)

    tke_sgs_ac = ps.AssignmentCollection(
        [lbm_method.conserved_quantity_computation.equilibrium_input_equations_from_pdfs(pdfs.center_vector),
         *tke_sgs_writer]
    )
    generate_sweep(ctx, "TurbulentChannel_TKE_SGS_Writer", tke_sgs_ac)

    #   Boundary conditions
    nu = lattice_viscosity_from_relaxation_rate(omega)
    u_tau_target = sp.Symbol("target_u_tau")

    noslip = NoSlip()
    freeslip_top = FreeSlip(stencil, normal_direction=normal_direction_top)

    wall_function_model = SpaldingsLaw(viscosity=nu, kappa=0.41, b=5.5, newton_steps=5)
    wfb_common_params = {
        'lb_method': lbm_method,
        'velocity_field': velocity,
        'wall_function_model': wall_function_model,
        'use_maronga_correction': False,
        'sampling_shift': 0,
        'filter_width': sp.Symbol("filter_width"),
        'data_type': data_type,
        'target_friction_velocity': u_tau_target
    }

    if wfb_weights_str == 'lattice':
        wfb_common_params['weight_method'] = WallFunctionBounce.WeightMethod.LATTICE_WEIGHT
    elif wfb_weights_str == 'geometric':
        wfb_common_params['weight_method'] = WallFunctionBounce.WeightMethod.GEOMETRIC_WEIGHT
    elif wfb_weights_str == 'none':
        pass
    else:
        raise ValueError("Invalid weight method")

    if velocity_str == 'mean':
        wfb_common_params['reference_velocity'] = WallFunctionBounce.ReferenceVelocity.MEAN_VELOCITY
    elif velocity_str == 'filtered':
        wfb_common_params['reference_velocity'] = WallFunctionBounce.ReferenceVelocity.FILTERED_VELOCITY
    else:
        raise ValueError("Invalid reference velocity")

    wfb_top = WallFunctionBounce(**wfb_common_params, normal_direction=normal_direction_top)
    wfb_bottom = WallFunctionBounce(**wfb_common_params, normal_direction=normal_direction_bottom)

    additional_data_handler_top = WFBAdditionalDataHandler(stencil, wfb_top, velocity, target=target)
    additional_data_handler_bottom = WFBAdditionalDataHandler(stencil, wfb_bottom, velocity, target=target)
    generate_boundary(ctx, "TurbulentChannel_NoSlip", noslip, lbm_method, target=target)
    generate_boundary(ctx, "TurbulentChannel_FreeSlip_top", freeslip_top, lbm_method, target=target)
    generate_boundary(ctx, "TurbulentChannel_WFB_bottom", wfb_bottom, lbm_method, target=target,
                      additional_data_handler=additional_data_handler_bottom)
    generate_boundary(ctx, "TurbulentChannel_WFB_top", wfb_top, lbm_method, target=target,
                      additional_data_handler=additional_data_handler_top)

    config_definitions = []
    if wfb_common_params['use_maronga_correction']:
        config_definitions.append("#define RUN_WITH_MARONGA_CORRECTION")
    if lbm_config.subgrid_scale_model is not None:
        config_definitions.append("#define RUN_WITH_SGS")
    if velocity_str == 'mean':
        config_definitions.append("#define MEAN_VELOCITY")
    elif velocity_str == 'filtered':
        config_definitions.append("#define FILTERED_VELOCITY")

    info_header_params = {
        'layout': layout,
        'd': stencil.D,
        'q': stencil.Q,
        'flow_axis': flow_axis,
        'wall_axis': wall_axis,
        'config_definitions': "\n".join(config_definitions)
    }

    ctx.write_file("CodegenIncludes.h", info_header.format(**info_header_params))
