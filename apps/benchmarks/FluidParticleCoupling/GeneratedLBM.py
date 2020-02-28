import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method_from_existing, create_lb_method
from lbmpy_walberla import generate_lattice_model
from pystencils_walberla import CodeGeneration

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.moments import is_even, get_order, MOMENT_SYMBOLS
from lbmpy.stencils import get_stencil


with CodeGeneration() as ctx:
    stencil = get_stencil("D3Q19", 'walberla')
    omega = sp.symbols("omega_:%d" % len(stencil))

    x, y, z = MOMENT_SYMBOLS
    one = sp.Rational(1, 1)
    sq = x ** 2 + y ** 2 + z ** 2
    moments = [
        [one, x, y, z],  # [0, 3, 5, 7]
        [sq - 1],  # [1]
        [3 * sq ** 2 - 6 * sq + 1],  # [2]
        [(3 * sq - 5) * x, (3 * sq - 5) * y, (3 * sq - 5) * z],  # [4, 6, 8]
        [3 * x ** 2 - sq, y ** 2 - z ** 2, x * y, y * z, x * z],  # [9, 11, 13, 14, 15]
        [(2 * sq - 3) * (3 * x ** 2 - sq), (2 * sq - 3) * (y ** 2 - z ** 2)],  # [10, 12]
        [(y ** 2 - z ** 2) * x, (z ** 2 - x ** 2) * y, (x ** 2 - y ** 2) * z]  # [16, 17, 18]
    ]
    method = create_lb_method(stencil=stencil, method='mrt', maxwellian_moments=False, nested_moments=moments)

    def modification_func(moment, eq, rate):
        omegaVisc = sp.Symbol("omega_visc")
        omegaBulk = ps.fields("omega_bulk: [3D]", layout='fzyx')
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




