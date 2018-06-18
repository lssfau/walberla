import sympy as sp
from lbmpy.boundaries import NoSlip, UBB
from lbmpy_walberla import Field, generate_lattice_model_files, RefinementScaling
from lbmpy.creationfunctions import create_lb_method
from lbmpy_walberla.boundary import create_boundary_class
from pystencils_walberla.cmake_integration import codegen
import pystencils as ps

# ------------- Lattice Model ------------------------------
force_field = ps.fields("force(3): [3D]", layout='fzyx')

omega = sp.Symbol("omega")

scaling = RefinementScaling()
scaling.add_standard_relaxation_rate_scaling(omega)
scaling.add_force_scaling(force_field)

generate_lattice_model_files(class_name='SrtWithForceFieldModel',
                             method='srt', stencil='D3Q19', force_model='guo', force=force_field.center_vector,
                             relaxation_rates=[omega], refinement_scaling=scaling)


def genBoundary():
    boundary = UBB([0.05, 0, 0], dim=3, name="MyUBB")
    method = create_lb_method(stencil='D3Q19', method='srt')
    return create_boundary_class(boundary, method)

def genNoSlip():
    boundary = NoSlip(name='MyNoSlip')
    method = create_lb_method(stencil='D3Q19', method='srt')
    return create_boundary_class(boundary, method)

codegen.register(['MyUBB.h', 'MyUBB.cpp'], genBoundary)
codegen.register(['MyNoSlip.h', 'MyNoSlip.cpp',], genNoSlip)

