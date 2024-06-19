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

    flow_axis = 0
    wall_axis = 1

    check_axis(flow_axis=flow_axis, wall_axis=wall_axis)

    #   ========================
    #      General Parameters
    #   ========================
    target = ps.Target.CPU

    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q19)
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
    sum_of_products = ps.fields(f"sum_of_products({stencil.D**2}): {data_type}[{stencil.D}D]", layout=layout)

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
                           subgrid_scale_model=SubgridScaleModel.QR,
                           # galilean_correction=True,
                           compressible=True,
                           omega_output_field=omega_field,
                           output={'velocity': velocity})

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
    # welford_update = welford_assignments(vector_field=velocity, mean_vector_field=mean_velocity)
    welford_update = welford_assignments(field=velocity, mean_field=mean_velocity,
                                         sum_of_products_field=sum_of_products)
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
    wfb_top = WallFunctionBounce(lbm_method, pdfs, normal_direction=normal_direction_top,
                                 wall_function_model=SpaldingsLaw(viscosity=nu,
                                                                  kappa=0.41, b=5.5, newton_steps=5),
                                 mean_velocity=mean_velocity, data_type=data_type,
                                 target_friction_velocity=u_tau_target)
    wfb_bottom = WallFunctionBounce(lbm_method, pdfs, normal_direction=normal_direction_bottom,
                                    wall_function_model=SpaldingsLaw(viscosity=nu,
                                                                     kappa=0.41, b=5.5, newton_steps=5),
                                    mean_velocity=mean_velocity, data_type=data_type,
                                    target_friction_velocity=u_tau_target)

    generate_boundary(ctx, "TurbulentChannel_NoSlip", noslip, lbm_method, target=target)
    generate_boundary(ctx, "TurbulentChannel_FreeSlip_top", freeslip_top, lbm_method, target=target)
    generate_boundary(ctx, "TurbulentChannel_WFB_bottom", wfb_bottom, lbm_method, target=target)
    generate_boundary(ctx, "TurbulentChannel_WFB_top", wfb_top, lbm_method, target=target)

    info_header_params = {
        'layout': layout,
        'd': stencil.D,
        'q': stencil.Q,
        'flow_axis': flow_axis,
        'wall_axis': wall_axis
    }

    ctx.write_file("CodegenIncludes.h", info_header.format(**info_header_params))
