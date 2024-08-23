import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration,generate_info_header, generate_pack_info_from_kernel,generate_sweep
from lbmpy import  LBStencil, Stencil, LBMOptimisation, LBMConfig, Method
from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule, create_lb_method
from lbmpy_walberla import generate_lbm_package, lbm_boundary_generator, generate_lattice_model
from lbmpy.boundaries import NeumannByCopy
from lbmpy.macroscopic_value_kernels import macroscopic_values_getter
from lbmpy import enums

info_header =   """
                const char * infoStencil = "{stencil}";
                const char * infoCollisionSetup = "{collision_setup}";
                const bool infoCseGlobal = {cse_global};
                const bool infoCsePdfs = {cse_pdfs};
                """

def get_relaxation_rates(stencil, omega):
    if stencil.Q == 19:
        return [omega] * 9 + [1.0] * 10
    elif stencil.Q == 27:
        return [omega] * 9 + [1.0] * 18
    
    return []

with CodeGeneration() as ctx:
    data_type   = "float64" if ctx.double_accuracy else "float32"
    stencil     = LBStencil(Stencil.D3Q27)
    method      = Method.MONOMIAL_CUMULANT

    omega       = sp.Symbol('omega')

    target      = ps.Target.CPU

    layout = 'fzyx'

    #   PDF Fields
    pdfs, pdfs_tmp = ps.fields(f'pdfs({stencil.Q}), pdfs_tmp({stencil.Q}): {data_type}[{stencil.D}D]', layout=layout)

    #   Velocity Field @ t = t+1
    u = ps.fields(f"u({stencil.D}): {data_type}[{stencil.D}D]", layout=layout)

    #   Density (scalar) Output Field
    c = ps.fields(f"c: {data_type}[{stencil.D}D]", layout=layout)

    macroscopic_fields = {  'density': c  }

    # LBM Optimisation
    lbm_opt = LBMOptimisation(cse_global=True,
                              symbolic_field=pdfs,
                              symbolic_temporary_field=pdfs_tmp,
                              field_layout=layout)

    #   ==================
    #      Method Setup
    #   ==================
    #! NOTE: Must set ω0 = ω1 for cumulant method or else velocity contribution to polynomial cumulants is zero
    #        causing an errenuous solution

    is_cumulant_method = (method == Method.SRT)

    relaxation_rates = get_relaxation_rates(stencil, omega) if is_cumulant_method else None

    lbm_config = LBMConfig(
        stencil=stencil,
        method=method,
        relaxation_rates=relaxation_rates if is_cumulant_method else None,
        relaxation_rate=omega if not is_cumulant_method else None,
        compressible=True,
        #streaming_pattern='esotwist',
        velocity_input=u
    )
    # subgrid_scale_model=enums.SubgridScaleModel.QR

    lb_method = create_lb_method(lbm_config=lbm_config)
    lbm_update_rule = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
   
    # Scale calculated lattice velocity to physical velocity for FDM routine
    print(collision_rule)
    
    generate_pack_info_from_kernel(ctx, "PackInfo", lbm_update_rule, target=target)

    neumann = lbm_boundary_generator(class_name='Neumann', flag_uid='Neumann',
                                    boundary_object=NeumannByCopy(),
                                    field_data_type=data_type)

    generate_lbm_package(ctx, name="AdvectionDiffusion",
                         collision_rule=collision_rule,
                         lbm_config=lbm_config, lbm_optimisation=lbm_opt, nonuniform=False,
                         boundaries=[neumann],
                         macroscopic_fields=macroscopic_fields,
                         data_type=data_type, pdfs_data_type=data_type)

    strategy = ps.simp.SimplificationStrategy()
    strategy.add(ps.simp.apply_to_all_assignments(sp.expand))
    strategy.add(ps.simp.sympy_cse)

    print(strategy.create_simplification_report(collision_rule))
    

    scale_factor = sp.Symbol("scale_factor")
    @ps.kernel
    def convertToSI():
        c.center @= c.center + scale_factor
    
    convertLBMtoSI = ps.AssignmentCollection( convertToSI, subexpressions=[] )

    generate_sweep(ctx, "SIConverter", convertLBMtoSI, target=target)

    infoHeaderParams = {
        'stencil': f"D{stencil.D}Q{stencil.Q}",
        'collision_setup': "CUMULANT",
        'cse_global': int(lbm_opt.cse_global),
        'cse_pdfs': int(lbm_opt.cse_pdfs),
    }
    
        
    generate_info_header(ctx, 'InfoHeader',
                         additional_code=info_header.format(**infoHeaderParams))
    