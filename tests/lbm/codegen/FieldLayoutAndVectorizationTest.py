import sympy as sp
from lbmpy.creationfunctions import create_lb_collision_rule
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model
from pystencils_walberla import get_vectorize_instruction_set

from collections import namedtuple

with CodeGeneration() as ctx:
    omega_shear = sp.symbols("omega")
    collision_rule = create_lb_collision_rule(stencil='D2Q9', compressible=False, method='srt')

    SetupDefinition = namedtuple('SetupDefinition', ['name', 'field_layout', 'vectorization_dict'])

    default_vectorize_instruction_set = get_vectorize_instruction_set(ctx)

    configurations = [SetupDefinition('FZYX_Vec', 'fzyx', {'instruction_set': default_vectorize_instruction_set}),
                      SetupDefinition('FZYX_NoVec', 'fzyx', {'instruction_set': None}),
                      SetupDefinition('ZYXF_Vec', 'zyxf', {'instruction_set': default_vectorize_instruction_set}),
                      # does/should not vectorize, but instead yield warning
                      SetupDefinition('ZYXF_NoVec', 'zyxf', {'instruction_set': None})]

    for conf in configurations:
        generate_lattice_model(ctx, 'FieldLayoutAndVectorizationTest_' + conf.name + '_LatticeModel', collision_rule,
                               field_layout=conf.field_layout, refinement_scaling=None,
                               cpu_vectorize_info=conf.vectorization_dict)
