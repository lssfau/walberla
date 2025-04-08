import sympy as sp
import numpy as np
import pystencils as ps

from lbmpy.enums import SubgridScaleModel
from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil

from lbmpy.creationfunctions import create_lb_update_rule, create_lb_collision_rule, create_lb_method

from lbmpy.boundaries import NoSlip, FreeSlip, UBB, ExtrapolationOutflow, NoSlipLinearBouzidi, QuadraticBounceBack

from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header
from lbmpy_walberla import generate_boundary, generate_lattice_model, generate_lbm_package, lbm_boundary_generator
from lbmpy_walberla.additional_data_handler import QuadraticBounceBackAdditionalDataHandler
from pystencils.typing import TypedSymbol

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q27)
    omega = sp.Symbol('omega')

    if ctx.gpu:
        target = ps.Target.GPU
    else:
        target = ps.Target.CPU

    layout = 'fzyx'

    compile_time_block_size = False

    max_threads = 256 if target == ps.Target.GPU else None

    if compile_time_block_size:
        sweep_block_size = (128, 1, 1)
    else:
        sweep_block_size = (TypedSymbol("cudaBlockSize0", np.int32),
                            TypedSymbol("cudaBlockSize1", np.int32),
                            TypedSymbol("cudaBlockSize2", np.int32))

    gpu_indexing_params = {'block_size': sweep_block_size}
    
    #   PDF Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}),pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]', layout=layout)

    #  Output Fields
    velocity = ps.fields(f"velocity({stencil.D}):{data_type}[{stencil.D}D]", layout=layout)
    density  = ps.fields(f"density({1}):{data_type}[{stencil.D}D]", layout=layout)

    macroscopic_fields = {'density': density, 'velocity': velocity}

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
                           galilean_correction=True,
                           zero_centered=True,
                           relaxation_rate=omega,
                           fourth_order_correction=0.01,
                           compressible=True,
                           output=macroscopic_fields)


    lbm_method = create_lb_method(lbm_config=lbm_config)
    lbm_update_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    collision_rule  = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    freeslip = lbm_boundary_generator("FreeSlipBC", flag_uid="FreeSlip", boundary_object=FreeSlip(stencil))

    noslip = lbm_boundary_generator(class_name='NoSlipBC', flag_uid='NoSlip', boundary_object=NoSlip(),
                                    field_data_type=data_type)
    
    noslipQbb = lbm_boundary_generator(class_name='NoSlipQBBBC', flag_uid='NoSlipQBB',
                                    boundary_object=QuadraticBounceBack( relaxation_rate=omega ),
                                    field_data_type=data_type)
    
    wallDistanceCallback = lambda *args: None
    obj_noslipQbb = lbm_boundary_generator(class_name='ObjNoSlipQBBBC', flag_uid='ObjNoSlipQBB',
                                    boundary_object=QuadraticBounceBack( relaxation_rate=omega ),
                                    field_data_type=data_type)

    outflow = lbm_boundary_generator( class_name='OutflowBC', flag_uid='Outflow',
                                      boundary_object=ExtrapolationOutflow(stencil[4], lbm_method),
                                      field_data_type=data_type)

    ubb = lbm_boundary_generator(class_name='UBBBC', flag_uid='UBB',
                                 boundary_object=UBB(lambda *args: None, dim=stencil.D, data_type=data_type),
                                 field_data_type=data_type)

    generate_lbm_package(ctx, name="LBComplexGeometry",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=False, boundaries=[outflow, ubb, freeslip, noslip, noslipQbb, obj_noslipQbb],
                         macroscopic_fields=macroscopic_fields,
                         gpu_indexing_params=gpu_indexing_params,
                         max_threads=max_threads, target=target,
                         data_type=data_type, pdfs_data_type=data_type,
                         cpu_openmp=ctx.openmp, set_pre_collision_pdfs=False)

    field_typedefs = {'VelocityField_T': velocity,
                      'ScalarField_T': density}

    stencil_typedefs = {'Stencil_T': stencil}

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'InfoHeader',stencil_typedefs=stencil_typedefs,
                         field_typedefs=field_typedefs)
