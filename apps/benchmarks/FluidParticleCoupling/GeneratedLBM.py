import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method_from_existing, create_lb_ast, create_lb_method
from lbmpy_walberla import generate_lattice_model
from pystencils_walberla import CodeGeneration

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.moments import MOMENT_SYMBOLS, is_even, get_order
from lbmpy.stencils import get_stencil
from lbmpy.maxwellian_equilibrium import get_moments_of_discrete_maxwellian_equilibrium
from lbmpy.forcemodels import *

from collections import OrderedDict
from lbmpy.methods import create_with_discrete_maxwellian_eq_moments

with CodeGeneration() as ctx:

    # methods: TRTlike, KBC, SRT, cumulant
    generatedMethod = "TRTlike"

    print("Generating " + generatedMethod + " LBM method")

    #
    # x,y,z = MOMENT_SYMBOLS
    #
    # e_sq = x**2 + y**2 + z**2
    #
    # moments = [0] * 19
    #
    # moments[0] = sp.sympify(1)
    # moments[1] = 19 * e_sq - 30
    # moments[2] = (21 * e_sq**2 - 53 * e_sq + 24) / 2
    #
    # moments[3] = x
    # moments[5] = y
    # moments[7] = z
    #
    # moments[4] = (5*e_sq - 9) * x
    # moments[6] = (5*e_sq - 9) * y
    # moments[8] = (5*e_sq - 9) * z
    #
    # moments[9] = 3 * x**2 - e_sq
    # moments[11] = y**2 - z**2
    #
    # moments[13] = x * y
    # moments[14] = y * z
    # moments[15] = x * z
    #
    # moments[10] = (3* e_sq - 5)*(3*x**2 - e_sq)
    # moments[12] = (3* e_sq - 5)*(y**2 - z**2)
    # moments[16] = (y**2 - z**2) * x
    # moments[17] = (z**2 - x**2) * y
    # moments[18] = (x**2 - y**2) * z
    #
    # m_eq = get_moments_of_discrete_maxwellian_equilibrium(stencil, moments, order=2, c_s_sq=sp.Rational(1,3))
    #
    # rr_dict = OrderedDict(zip(moments, omega))
    #
    # forcing = sp.symbols("forcing_:%d" % 3)
    # forcing=(sp.symbols("fx"),0,0)
    # forcemodel=Luo(forcing) #None
    # method = create_with_discrete_maxwellian_eq_moments(stencil, rr_dict, compressible=False, force_model=forcemodel)
    # at this stage, we have the MRT model from the LBM book!

     #this is the schiller / walberla MRT

    if generatedMethod == "TRTlike":
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

         # optimizations

         collision_rule = create_lb_collision_rule(lb_method=my_method,
                             optimization={"cse_global": True,
                                           "cse_pdfs": False,
                                           #"split": True
                                           } )

         generate_lattice_model(ctx, 'GeneratedLBM', collision_rule)

    elif generatedMethod == "KBC":
        methodName = 'trt-kbc-n4'
        #omega_field = ps.fields("omegaField(1): [3D]", layout='fzyx')  omega_output_field=omega_field.center,
        collision_rule = create_lb_collision_rule(method=methodName,entropic=True, stencil=get_stencil("D3Q27", 'walberla'), compressible=True,
                        optimization={"cse_global": True,
                                      "cse_pdfs": False,
                                      "split": True})
        generate_lattice_model(ctx, 'GeneratedLBM', collision_rule)

    elif generatedMethod == "SRT":
        methodName = 'srt'
        collision_rule = create_lb_collision_rule(method=methodName,stencil=get_stencil("D3Q19", 'walberla'),
                     optimization={"cse_global": True,
                                   "cse_pdfs": False,
                                   "split": True})
        generate_lattice_model(ctx, 'GeneratedLBM', collision_rule)

    elif generatedMethod == "cumulant":

        x,y,z = MOMENT_SYMBOLS

        cumulants = [0] * 27

        cumulants[0] = sp.sympify(1) #000

        cumulants[1] = x #100
        cumulants[2] = y #010
        cumulants[3] = z #001

        cumulants[4] = x*y #110
        cumulants[5] = x*z #101
        cumulants[6] = y*z #011

        cumulants[7] = x**2 - y**2 #200 - 020
        cumulants[8] = x**2 - z**2 #200 - 002
        cumulants[9] = x**2 + y**2 + z**2 #200 + 020 + 002

        cumulants[10] = x*y**2 + x*z**2 #120 + 102
        cumulants[11] = x**2*y + y*z**2 #210 + 012
        cumulants[12] = x**2*z + y**2*z #201 + 021
        cumulants[13] = x*y**2 - x*z**2 #120 - 102
        cumulants[14] = x**2*y - y*z**2 #210 - 012
        cumulants[15] = x**2*z - y**2*z #201 - 021

        cumulants[16] = x*y*z #111

        cumulants[17] = x**2*y**2 - 2*x**2*z**2 + y**2*z**2 # 220- 2*202 +022
        cumulants[18] = x**2*y**2 + x**2*z**2 - 2*y**2*z**2 # 220 + 202 - 2*002
        cumulants[19] = x**2*y**2 + x**2*z**2 + y**2*z**2 # 220 + 202 + 022

        cumulants[20] = x**2*y*z # 211
        cumulants[21] = x*y**2*z # 121
        cumulants[22] = x*y*z**2 # 112

        cumulants[23] = x**2*y**2*z # 221
        cumulants[24] = x**2*y*z**2 # 212
        cumulants[25] = x*y**2*z**2 # 122

        cumulants[26] = x**2*y**2*z**2 # 222


        from lbmpy.moments import get_order
        def get_relaxation_rate(cumulant, omega):
            if get_order(cumulant) <= 1:
                return 0
            elif get_order(cumulant) == 2 and cumulant != x**2 + y**2 + z**2:
                return omega
            else:
                return 1

        stencil = get_stencil("D3Q27", 'walberla')

        omega = sp.Symbol("omega")
        rr_dict = OrderedDict((c, get_relaxation_rate(c, omega))
                              for c in cumulants)

        from lbmpy.methods import create_with_continuous_maxwellian_eq_moments
        my_method = create_with_continuous_maxwellian_eq_moments(stencil, rr_dict, cumulant=True, compressible=True)

        collision_rule = create_lb_collision_rule(lb_method=my_method,
                            optimization={"cse_global": True,
                                          "cse_pdfs": False})

        generate_lattice_model(ctx, 'GeneratedLBM', collision_rule)

    else:
        print("Invalid generated method string! " + generatedMethod)



