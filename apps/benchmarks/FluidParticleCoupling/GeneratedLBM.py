import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method_from_existing, create_lb_ast, create_lb_method
from lbmpy_walberla import generate_lattice_model
from pystencils_walberla import CodeGeneration

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.moments import MOMENT_SYMBOLS, is_even, get_order
from lbmpy.stencils import get_stencil
from lbmpy.forcemodels import *


with CodeGeneration() as ctx:
    stencil = get_stencil("D3Q19", 'walberla')
    omega = sp.symbols("omega_:%d" % len(stencil))
    method = create_lb_method(stencil=stencil, method='mrt', maxwellian_moments=False)

    def modification_func(moment, eq, rate):
        omegaVisc = sp.Symbol("omega_visc")
        omegaBulk = ps.fields("omega_bulk: [3D]", layout='fzyx')# sp.Symbol("omega_bulk")
        omegaMagic = sp.Symbol("omega_magic")
        if get_order(moment) <= 1:
            return moment, eq, 0
        elif rate == omega[1]:
            return moment, eq, omegaBulk.center_vector
        elif rate == omega[2]:
            return moment, eq, omegaBulk.center_vector
        elif is_even(moment):
            return moment, eq, omegaVisc
        else:
            return moment, eq, omegaMagic

    my_method = create_lb_method_from_existing(method, modification_func)

    collision_rule = create_lb_collision_rule(lb_method=my_method,
                     optimization={"cse_global": True,
                                   "cse_pdfs": False
                                   } )

    generate_lattice_model(ctx, 'GeneratedLBM', collision_rule)




