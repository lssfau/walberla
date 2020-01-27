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

    forcing=(sp.symbols("fx"),0,0)
    forcemodel=Luo(forcing) #None

    # methods: TRTlike, KBC, SRT
    generatedMethod = "TRTlike"

    if generatedMethod == "TRTlike":
        stencil = get_stencil("D3Q19", 'walberla')
        omega = sp.symbols("omega_:%d" % len(stencil))

        methodWithForce = create_lb_method(stencil=stencil, method='mrt', maxwellian_moments=False,force_model=forcemodel)

        def modification_func(moment, eq, rate):
            omegaVisc = sp.Symbol("omega_visc")
            omegaBulk = ps.fields("omega_bulk: [3D]", layout='fzyx')# = sp.Symbol("omega_bulk")
            #omegaBulk = sp.Symbol("omega_bulk")
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

        my_methodWithForce = create_lb_method_from_existing(methodWithForce, modification_func)

        collision_rule = create_lb_collision_rule(lb_method=my_methodWithForce)
        #    ,
        #                     optimization={"cse_global": True,
        #                                   "cse_pdfs": False,
        #                                   "split": True,
        #                                   "vectorization":{'instruction_set': 'avx',
        #                                                    'nontemporal': True,
        #                                                    'assume_inner_stride_one': True}})

        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule)

    elif generatedMethod == "SRT":
        collision_rule = create_lb_collision_rule(method='srt',stencil=get_stencil("D3Q19", 'walberla'), force_model=forcemodel,
                            optimization={"cse_global": True,
                                          "cse_pdfs": False,
                                          #"split": True
                                          })
        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule)

    elif generatedMethod == "TRT":
        collision_rule = create_lb_collision_rule(method='trt',stencil=get_stencil("D3Q19", 'walberla'), force_model=forcemodel, maxwellian_moments=False)
        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule)


    elif generatedMethod == "KBC":
        collision_rule = create_lb_collision_rule(method='trt-kbc-n4',entropic=True, stencil=get_stencil("D3Q27", 'walberla'), compressible=True, force_model=forcemodel,
                            optimization={"cse_global": True,
                                          "cse_pdfs": False,
                                          #"split": True
                                          })
        generate_lattice_model(ctx, 'GeneratedLBMWithForce', collision_rule)

    else:
        print("Invalid generated method string! " + generatedMethod)
