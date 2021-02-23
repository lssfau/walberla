import sympy as sp
from sympy.core.cache import clear_cache
import pystencils as ps
from lbmpy.creationfunctions import create_lb_method_from_existing, create_lb_ast, create_lb_method
from lbmpy_walberla import generate_lattice_model
from pystencils_walberla import CodeGeneration
from pystencils_walberla import get_vectorize_instruction_set

from lbmpy.creationfunctions import create_lb_collision_rule
from lbmpy.moments import MOMENT_SYMBOLS, is_even, get_order
from lbmpy.stencils import get_stencil
from lbmpy.forcemodels import Luo

from collections import OrderedDict

with CodeGeneration() as ctx:

    forcing=(sp.symbols("fx"),0,0)
    forcemodel=Luo(forcing)

    generatedMethod = "TRTlike"
    #generatedMethod = "D3Q27TRTlike"
    #generatedMethod = "cumulant"
    #generatedMethod = "cumulantTRT"

    print("Generating " + generatedMethod + " LBM method")

    clear_cache()

    cpu_vectorize_info = {'instruction_set': get_vectorize_instruction_set(ctx)}

    if generatedMethod == "TRTlike":
        omegaVisc = sp.Symbol("omega_visc")
        omegaBulk = ps.fields("omega_bulk: [3D]", layout='fzyx')
        omegaMagic = sp.Symbol("omega_magic")
        stencil = get_stencil("D3Q19", 'walberla')

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

        # relaxation rate for first group of moments (1,x,y,z) is set to zero internally
        relaxation_rates=[omegaBulk.center_vector, omegaBulk.center_vector, omegaMagic, omegaVisc, omegaVisc, omegaMagic]

        methodWithForce = create_lb_method(stencil=stencil, method='mrt', maxwellian_moments=False,
                                           force_model=forcemodel, nested_moments=moments, relaxation_rates=relaxation_rates)

        #print(methodWithForce.relaxation_rates)
        #print(methodWithForce.moment_matrix)
        collision_rule = create_lb_collision_rule(lb_method=methodWithForce, optimization={'cse_global': True})
        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule, field_layout='fzyx', cpu_vectorize_info=cpu_vectorize_info)
 
    if generatedMethod == "D3Q27TRTlike":

        omegaVisc = sp.Symbol("omega_visc")
        omegaBulk = ps.fields("omega_bulk: [3D]", layout='fzyx')
        omegaMagic = sp.Symbol("omega_magic")
        stencil = get_stencil("D3Q27", 'walberla')

        relaxation_rates=[omegaVisc, omegaBulk.center_vector, omegaMagic, omegaVisc, omegaMagic, omegaVisc]

        methodWithForce = create_lb_method(stencil=stencil, method='mrt', maxwellian_moments=False,weighted=True, compressible=False,
                                           force_model=forcemodel, relaxation_rates=relaxation_rates)

        collision_rule = create_lb_collision_rule(lb_method=methodWithForce, optimization={'cse_global': True})
        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule, field_layout='fzyx', cpu_vectorize_info=cpu_vectorize_info)

    if generatedMethod == "cumulant":

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
        my_method = create_with_continuous_maxwellian_eq_moments(stencil, rr_dict, cumulant=True, compressible=True, force_model=forcemodel)

        collision_rule = create_lb_collision_rule(lb_method=my_method,
                                                  optimization={"cse_global": True,
                                                                "cse_pdfs": False})

        print(my_method.relaxation_rates)

        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule, field_layout='fzyx', cpu_vectorize_info=cpu_vectorize_info)


    if generatedMethod == "cumulantTRT":

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

        def get_relaxation_rate(cumulant, omegaVisc, omegaMagic):
            if get_order(cumulant) <= 1:
                return 0
            elif is_even(cumulant):
                return omegaVisc
            else:
                return omegaMagic

        stencil = get_stencil("D3Q27", 'walberla')

        omegaVisc = sp.Symbol("omega_visc")
        omegaMagic = sp.Symbol("omega_magic")

        rr_dict = OrderedDict((c, get_relaxation_rate(c, omegaVisc, omegaMagic))
                              for c in cumulants)

        from lbmpy.methods import create_with_continuous_maxwellian_eq_moments
        my_method = create_with_continuous_maxwellian_eq_moments(stencil, rr_dict, cumulant=True, compressible=True, force_model=forcemodel)

        collision_rule = create_lb_collision_rule(lb_method=my_method,
                                                  optimization={"cse_global": True,
                                                                "cse_pdfs": False})

        print(my_method.relaxation_rates)

        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule, field_layout='fzyx', cpu_vectorize_info=cpu_vectorize_info)

