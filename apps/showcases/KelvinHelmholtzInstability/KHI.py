import sympy as sp
import pystencils as ps
import numpy as np

from lbmpy.enums import SubgridScaleModel
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.flow_statistics import welford_assignments
from lbmpy.relaxationrates import lattice_viscosity_from_relaxation_rate

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule

from lbmpy.boundaries import FreeSlip, NeumannByCopy

from pystencils import Assignment

from pystencils_walberla import CodeGeneration, generate_sweep_collection, generate_info_header

from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator

from pystencils.typing import TypedSymbol

info_header = """
namespace walberla {{

    const FlagUID FluidFlagUID("Fluid");
    
    using flag_t = uint16_t;
    using FlagField_T = FlagField< flag_t >;

    using StorageSpecification_T        = lbm::{lbm_name}StorageSpecification;
    using ScalarStorageSpecification_T  = lbm::{scalar_name}StorageSpecification;

    using LBMCommunicationStencil_T        = StorageSpecification_T::CommunicationStencil;
    using ScalarLBMCommunicationStencil_T  = ScalarStorageSpecification_T::CommunicationStencil;

    using PdfField_T                 = lbm_generated::PdfField< StorageSpecification_T >;
    using ScalarPdfField_T           = lbm_generated::PdfField< ScalarStorageSpecification_T >;

    using SweepCollection_T          = lbm::{lbm_name}SweepCollection;
    using ScalarSweepCollection_T    = lbm::{scalar_name}SweepCollection;

    #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
        using GPUField_T              = gpu::GPUField< real_t >;
        
        using GPUPdfField_T           = lbm_generated::GPUPdfField< StorageSpecification_T >;
        using PackInfo_T              = lbm_generated::UniformGeneratedGPUPdfPackInfo< GPUPdfField_T >;

        using ScalarGPUPdfField_T     = lbm_generated::GPUPdfField< ScalarStorageSpecification_T >;
        using ScalarPackInfo_T        = lbm_generated::UniformGeneratedGPUPdfPackInfo< ScalarGPUPdfField_T >;

        using gpu::communication::UniformGPUScheme;
        using Timeloop_T = DeviceSynchronizeSweepTimeloop;

    #else
        using PackInfo_T       = lbm_generated::UniformGeneratedPdfPackInfo< PdfField_T >;
        using ScalarPackInfo_T = lbm_generated::UniformGeneratedPdfPackInfo< ScalarPdfField_T >;
        
        using blockforest::communication::UniformBufferedScheme;
        using Timeloop_T = SweepTimeloop;
    #endif

    using BoundaryCollection_T = lbm::{lbm_name}BoundaryCollection< FlagField_T >;
    using ScalarBoundaryCollection_T = lbm::{scalar_name}BoundaryCollection< FlagField_T >;


    namespace codegen {{
        static const std::string configID = "{configID}";
        static constexpr field::Layout layout = field::{layout};
        static constexpr uint_t flowAxis = {flow_axis};
        static constexpr uint_t wallAxis = {wall_axis};
    }}
}}
"""

options_dict = {
    'smagorinsky': {
        'subgrid_scale_model': SubgridScaleModel.SMAGORINSKY,
        'fourth_order_correction': True,
    },
    'qr': {
        'subgrid_scale_model': SubgridScaleModel.QR,
        'fourth_order_correction': True,
    },
    'cumulant-allOne': {
        'fourth_order_correction': False,
    },
    'cumulant': {
        'fourth_order_correction': True,
    },
    'cumulant-L1p000': {
        'fourth_order_correction': 1.0
    },
    'cumulant-L0p100': {
        'fourth_order_correction': 0.1
    },
    'cumulant-L0p010': {
        'fourth_order_correction': 0.01
    },
    'cumulant-L0p001': {
        'fourth_order_correction': 0.001
    }
}


def check_axis(flow_axis, wall_axis):
    assert flow_axis != wall_axis, "Axes must be distinct."
    assert all(0 <= axis < 3 for axis in (flow_axis, wall_axis)), "Axes must be between 0 and 2."

with CodeGeneration() as ctx:
    config_tokens = ctx.config.split('_')
    config_handle = config_tokens[0]     # Set by CMakeLists.txt, Options:
                                         # 'smagorinsky','qr','cumulant','cumulant-L1p000',
                                         # 'cumulant-L0p100','cumulant-L0p010','cumulant-L0p001',
                                         # 'cumulant-allOne'

    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil         = LBStencil(Stencil.D3Q27)
    scalar_stencil  = LBStencil(Stencil.D3Q19)

    omega = sp.Symbol('omega')
    scalar_omega = sp.Symbol('omega')

    layout = 'fzyx'

    compile_time_block_size = False

    flow_axis  = 0
    wall_axis  = 1

    check_axis(flow_axis=flow_axis, wall_axis=wall_axis)

    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU

    max_threads = 256 if target == ps.Target.GPU else None

    if compile_time_block_size:
        sweep_block_size = (128, 1, 1)
    else:
        sweep_block_size = (TypedSymbol("gpuBlockSize0", np.int32),
                            TypedSymbol("gpuBlockSize1", np.int32),
                            TypedSymbol("gpuBlockSize2", np.int32))

    gpu_indexing_params = {'block_size': sweep_block_size}

    #--------------------------------------------------------------------------------------------------------------
    # Velocity Config
    #--------------------------------------------------------------------------------------------------------------
    #   PDF Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]', layout=layout)

    #   Output Fields
    density = ps.fields(f"density({1}): {data_type}[{stencil.D}D]", layout=layout)
    velocity = ps.fields(f"velocity({stencil.D}): {data_type}[{stencil.D}D]", layout=layout)

    macroscopic_fields = {'velocity': velocity,
                          'density': density}

    lbm_opt = LBMOptimisation(cse_global=True,
                              symbolic_field=pdfs,
                              symbolic_temporary_field=pdfs_tmp,
                              field_layout=layout)
    
    options = options_dict[config_handle]

    lbm_config = LBMConfig(stencil=stencil,
                           method=Method.CUMULANT,
                           relaxation_rate=omega,
                           galilean_correction=True,
                           compressible=True,
                           zero_centered=True,
                           output=macroscopic_fields,
                           **options)
              
    collision_rule  = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    freeslip = lbm_boundary_generator( class_name='FreeSlip_BC', flag_uid='FreeSlip', 
                                      boundary_object=FreeSlip(stencil), 
                                      field_data_type=data_type )
    lbm_name = "KHI"
    generate_lbm_package(ctx, name=lbm_name, collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=False, boundaries=[freeslip],
                         macroscopic_fields=macroscopic_fields,
                         gpu_indexing_params=gpu_indexing_params,
                         max_threads=max_threads, target=target,
                         data_type=data_type, pdfs_data_type=data_type,
                         cpu_openmp=ctx.openmp, set_pre_collision_pdfs=False)

    #--------------------------------------------------------------------------------------------------------------
    # Advection-Diffusion Config
    #--------------------------------------------------------------------------------------------------------------
    scalar_pdfs, scalar_pdfs_tmp = ps.fields(f'pdfs({scalar_stencil.Q}), pdfs_tmp({scalar_stencil.Q}): {data_type}[{scalar_stencil.D}D]', layout=layout)
    scalar = ps.fields(f"scalar({1}): {data_type}[{scalar_stencil.D}D]", layout=layout)

    macroscopic_fields = {'density': scalar}

    scalar_lbm_opt = LBMOptimisation(cse_global=True,
                                    symbolic_field=scalar_pdfs,
                                    symbolic_temporary_field=scalar_pdfs_tmp,
                                    field_layout=layout)

    scalar_lbm_config = LBMConfig(  stencil=scalar_stencil,
                                    method=Method.MRT,
                                    relaxation_rate=omega,
                                    compressible=True,
                                    zero_centered=True,
                                    velocity_input=velocity,
                                    output=macroscopic_fields )

    scalar_update_rule     = create_lb_update_rule(lbm_config=scalar_lbm_config, lbm_optimisation=scalar_lbm_opt)
    scalar_collision_rule  = create_lb_collision_rule(lbm_config=scalar_lbm_config, lbm_optimisation=scalar_lbm_opt)
    scalar_lbm_method = scalar_update_rule.method

    neumann = lbm_boundary_generator( class_name='Neumann_BC', flag_uid='Neumann', 
                                      boundary_object=NeumannByCopy(scalar_lbm_method), 
                                      field_data_type=data_type )

    scalar_name = "KHI_AdvectionDiffusion"
    generate_lbm_package(ctx, name=scalar_name,collision_rule=scalar_collision_rule,
                         lbm_config=scalar_lbm_config, lbm_optimisation=scalar_lbm_opt,
                         nonuniform=False, boundaries=[neumann],
                         macroscopic_fields=macroscopic_fields,
                         gpu_indexing_params=gpu_indexing_params,
                         max_threads=max_threads, target=target,
                         data_type=data_type, pdfs_data_type=data_type,
                         cpu_openmp=ctx.openmp, set_pre_collision_pdfs=False)


    info_header_params = {
        'configID':config_handle,
        'lbm_name': lbm_name,
        'scalar_name':scalar_name,
        'layout': layout,
        'flow_axis': flow_axis,
        'wall_axis': wall_axis
    }

    field_typedefs = {'VectorField_T': velocity,
                      'ScalarField_T': density }

    stencil_typedefs = {'Stencil_T': stencil,
                        'ScalarStencil_T': scalar_stencil}

    generate_info_header(ctx, 'InfoHeader',stencil_typedefs=stencil_typedefs,
                         field_typedefs=field_typedefs, additional_headers={"field/FlagField.h"},
                         additional_code=info_header.format(**info_header_params))
