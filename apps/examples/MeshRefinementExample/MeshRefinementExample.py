import sympy as sp
import pystencils as ps

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_method

from lbmpy.boundaries import NoSlip, FreeSlip, UBB, ExtrapolationOutflow

from pystencils_walberla import CodeGeneration, generate_info_header
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator


with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    stencil = LBStencil(Stencil.D3Q19)
    omega = sp.Symbol('omega')

    target = ps.Target.CPU
    layout = 'fzyx'

    #   Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}),pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]', layout=layout)
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
    lbm_config = LBMConfig(stencil=stencil, output=macroscopic_fields)


    lbm_method = create_lb_method(lbm_config=lbm_config)
    collision_rule  = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)

    freeslip = lbm_boundary_generator("FreeSlipBC", flag_uid="FreeSlip", boundary_object=FreeSlip(stencil))

    noslip = lbm_boundary_generator(class_name='NoSlipBC', flag_uid='NoSlip',
                                    boundary_object=NoSlip(calculate_force_on_boundary=True),
                                    field_data_type=data_type)

    outflow = lbm_boundary_generator( class_name='OutflowBC', flag_uid='Outflow',
                                      boundary_object=ExtrapolationOutflow(stencil[4], lbm_method),
                                      field_data_type=data_type)

    inlet_velocity = (sp.symbols("u_x"), 0, 0) if stencil.D == 3 else (sp.symbols("u_x"), 0)
    ubb = lbm_boundary_generator(class_name='UBBBC', flag_uid='UBB',
                                 boundary_object=UBB(inlet_velocity, density=1.0, data_type=data_type, dim=stencil.D),
                                 field_data_type=data_type)

    generate_lbm_package(ctx, name="MeshRefinementExample",
                         collision_rule=collision_rule, lbm_config=lbm_config, lbm_optimisation=lbm_opt,
                         nonuniform=True, boundaries=[outflow, ubb, noslip, freeslip],
                         macroscopic_fields=macroscopic_fields, target=target, data_type=data_type,
                         pdfs_data_type=data_type, cpu_openmp=ctx.openmp)

    field_typedefs = {'VelocityField_T': velocity, 'ScalarField_T': density}
    stencil_typedefs = {'Stencil_T': stencil}
    generate_info_header(ctx, 'InfoHeader', stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
