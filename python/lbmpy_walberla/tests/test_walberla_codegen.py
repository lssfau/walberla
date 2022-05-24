import unittest

import sympy as sp

import pystencils as ps
from lbmpy import ForceModel, LBMConfig, LBMOptimisation, LBStencil, Method, Stencil
from lbmpy.boundaries import UBB, NoSlip
from lbmpy.creationfunctions import create_lb_method, create_lb_update_rule, create_lb_collision_rule
from lbmpy_walberla import RefinementScaling, generate_boundary, generate_lattice_model
from lbmpy_walberla.sparse import ListLbGenerator
from pystencils_walberla import generate_pack_info_for_field, generate_pack_info_from_kernel
from pystencils_walberla.cmake_integration import ManualCodeGenerationContext


class WalberlaLbmpyCodegenTest(unittest.TestCase):

    @staticmethod
    def test_lattice_model():
        with ManualCodeGenerationContext() as ctx:
            force_field = ps.fields("force(3): [3D]", layout='fzyx')
            omega = sp.Symbol("omega")

            stencil = LBStencil(Stencil.D3Q19)
            lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rates=[omega], compressible=True,
                                   force_model=ForceModel.GUO, force=force_field.center_vector)

            cr = create_lb_collision_rule(lbm_config=lbm_config)

            scaling = RefinementScaling()
            scaling.add_standard_relaxation_rate_scaling(omega)
            scaling.add_force_scaling(force_field)

            generate_lattice_model(ctx, 'SrtWithForceFieldModel', cr, refinement_scaling=scaling)
            generate_boundary(ctx, 'MyUBB', UBB([0.05, 0, 0]), cr.method)
            generate_boundary(ctx, 'MyNoSlip', NoSlip(), cr.method)
            assert 'static const bool compressible = true;' in ctx.files['SrtWithForceFieldModel.h']

    @staticmethod
    def test_sparse():
        from lbmpy.creationfunctions import create_lb_collision_rule
        from pystencils import get_code_str
        g = ListLbGenerator(create_lb_collision_rule())
        kernel_code = get_code_str(g.kernel())
        assert 'num_cells' in kernel_code
        setter_code = get_code_str(g.setter_ast())
        assert 'num_cells' in setter_code
        getter_code = get_code_str(g.getter_ast())
        assert 'num_cells' in getter_code

    @staticmethod
    def test_pack_info():
        with ManualCodeGenerationContext() as ctx:
            f = ps.fields("f(9): [3D]")
            generate_pack_info_for_field(ctx, 'MyPackInfo1', f)

            stencil = LBStencil(stencil=Stencil.D3Q19)
            lbm_config = LBMConfig(stencil=stencil, method=Method.SRT)

            lb_assignments = create_lb_update_rule(lbm_config=lbm_config).main_assignments
            generate_pack_info_from_kernel(ctx, 'MyPackInfo2', lb_assignments)

    @staticmethod
    def test_incompressible():
        with ManualCodeGenerationContext() as ctx:
            omega = sp.Symbol("omega")

            stencil = LBStencil(stencil=Stencil.D3Q19)
            lbm_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rates=[omega], compressible=False)

            cr = create_lb_collision_rule(lbm_config=lbm_config)
            generate_lattice_model(ctx, 'Model', cr)
            assert 'static const bool compressible = false;' in ctx.files['Model.h']

    @staticmethod
    def test_output_field():
        with ManualCodeGenerationContext(openmp=True, double_accuracy=True) as ctx:
            omega_field = ps.fields("omega_out: [3D]", layout='fzyx')

            stencil = LBStencil(stencil=Stencil.D3Q27)
            lbm_config = LBMConfig(stencil=stencil, method=Method.TRT_KBC_N1,
                                   zero_centered=False,
                                   relaxation_rates=[1.5, sp.Symbol('omega_free')],
                                   compressible=True, entropic=True,
                                   omega_output_field=omega_field)

            cr = create_lb_collision_rule(lbm_config=lbm_config)
            generate_lattice_model(ctx, 'Model', cr)

    @staticmethod
    def test_fluctuating():
        with ManualCodeGenerationContext(openmp=True, double_accuracy=True) as ctx:
            omega_shear = sp.symbols("omega")
            force_field, vel_field = ps.fields("force(3), velocity(3): [3D]", layout='fzyx')

            stencil = LBStencil(stencil=Stencil.D3Q19)
            lbm_config = LBMConfig(stencil=stencil, method=Method.MRT, compressible=True, zero_centered=False,
                                   fluctuating={'seed': 0, 'temperature': 1e-6}, force_model=ForceModel.GUO,
                                   relaxation_rates=[omega_shear] * 19, force=force_field.center_vector)

            lbm_opt = LBMOptimisation(cse_global=True)

            # the collision rule of the LB method where the some advanced features
            collision_rule = create_lb_collision_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
            generate_lattice_model(ctx, 'FluctuatingMRT', collision_rule)

    @staticmethod
    def test_boundary_3D():
        with ManualCodeGenerationContext(openmp=True, double_accuracy=True, cuda=True) as ctx:
            stencil = LBStencil(stencil=Stencil.D3Q19)
            lbm_config = LBMConfig(stencil=stencil, method=Method.SRT)

            lb_method = create_lb_method(lbm_config=lbm_config)
            generate_boundary(ctx, 'Boundary', NoSlip(), lb_method, target=ps.Target.GPU)

    @staticmethod
    def test_boundary_2D():
        with ManualCodeGenerationContext(openmp=True, double_accuracy=True, cuda=True) as ctx:
            stencil = LBStencil(stencil=Stencil.D2Q9)
            lbm_config = LBMConfig(stencil=stencil, method=Method.SRT)

            lb_method = create_lb_method(lbm_config=lbm_config)
            generate_boundary(ctx, 'Boundary', NoSlip(), lb_method, target=ps.Target.GPU)


if __name__ == '__main__':
    unittest.main()
