# import warnings
from typing import Callable, List


import numpy as np
import sympy as sp
from jinja2 import Environment, PackageLoader, StrictUndefined, Template
from sympy.tensor import IndexedBase

import pystencils as ps
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor, StreamPullTwoFieldsAccessor
from lbmpy.relaxationrates import relaxation_rate_scaling
from lbmpy.updatekernels import create_lbm_kernel, create_stream_only_kernel
from pystencils import AssignmentCollection, create_kernel, Target
from pystencils.astnodes import SympyAssignment
from pystencils.backends.cbackend import CBackend, CustomSympyPrinter, get_headers
from pystencils.typing import BasicType, CastFunc, TypedSymbol
from pystencils.field import Field
from pystencils.node_collection import NodeCollection
from pystencils.stencil import offset_to_direction_string
from pystencils.sympyextensions import get_symmetric_part
from pystencils.typing.transformations import add_types

from pystencils_walberla.kernel_info import KernelInfo
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.utility import config_from_context

cpp_printer = CustomSympyPrinter()
REFINEMENT_SCALE_FACTOR = sp.Symbol("level_scale_factor")


def __type_equilibrium_assignments(assignments, config, subs_dict):
    result = assignments.new_with_substitutions(subs_dict)
    result = NodeCollection(result.main_assignments)
    result.evaluate_terms()
    result = add_types(result.all_assignments, config)
    return result


def __lattice_model(generation_context, class_name, config, lb_method, stream_collide_ast, collide_ast, stream_ast,
                    refinement_scaling):
    stencil_name = lb_method.stencil.name
    if not stencil_name:
        raise ValueError("lb_method uses a stencil that is not supported in waLBerla")

    communication_stencil_name = stencil_name if stencil_name != "D3Q15" else "D3Q27"
    default_dtype = config.data_type.default_factory()

    cqc = lb_method.conserved_quantity_computation

    vel_symbols = cqc.velocity_symbols
    rho_sym = sp.Symbol('rho')
    reference_density = rho_sym if cqc.compressible else cqc.background_density
    pdfs_sym = sp.symbols(f'f_:{lb_method.stencil.Q}')

    vel_arr_symbols = [IndexedBase(TypedSymbol('u', default_dtype), shape=(1,))[i] for i in range(len(vel_symbols))]
    subs_dict = {a: b for a, b in zip(vel_symbols, vel_arr_symbols)}

    momentum_density_symbols = sp.symbols(f'md_:{len(vel_symbols)}')

    equilibrium = lb_method.get_equilibrium()
    lhs_list = [a.lhs for a in equilibrium.main_assignments]

    equilibrium_matrix = sp.Matrix([e.rhs for e in equilibrium.main_assignments])
    symmetric_equilibrium_matrix = get_symmetric_part(equilibrium_matrix, vel_symbols)
    asymmetric_equilibrium_matrix = sp.expand(equilibrium_matrix - symmetric_equilibrium_matrix)

    equilibrium = AssignmentCollection([ps.Assignment(lhs, rhs)
                                        for lhs, rhs in zip(lhs_list, equilibrium_matrix)])
    symmetric_equilibrium = AssignmentCollection([ps.Assignment(lhs, rhs)
                                                  for lhs, rhs in zip(lhs_list, symmetric_equilibrium_matrix)])
    asymmetric_equilibrium = AssignmentCollection([ps.Assignment(lhs, rhs)
                                                   for lhs, rhs in zip(lhs_list, asymmetric_equilibrium_matrix)])

    equilibrium = __type_equilibrium_assignments(equilibrium, config, subs_dict)
    symmetric_equilibrium = __type_equilibrium_assignments(symmetric_equilibrium, config, subs_dict)
    asymmetric_equilibrium = __type_equilibrium_assignments(asymmetric_equilibrium, config, subs_dict)

    force_model = lb_method.force_model
    macroscopic_velocity_shift = None
    if force_model:
        if hasattr(force_model, 'macroscopic_velocity_shift'):
            macroscopic_velocity_shift = [e.subs(force_model.subs_dict_force)
                                          for e in force_model.macroscopic_velocity_shift(reference_density)]
            macroscopic_velocity_shift = [expression_to_code(e, "lm.", ['rho'], dtype=default_dtype)
                                          for e in macroscopic_velocity_shift]

    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(sp.Symbol('rho_in'), vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(eq_input_from_input_eqs, dtype=default_dtype,
                                                                   variables_without_prefix=['rho_in', 'u'])
    momentum_density_getter = cqc.output_equations_from_pdfs(pdfs_sym, {'density': rho_sym,
                                                                        'momentum_density': momentum_density_symbols})

    is_float = True if issubclass(default_dtype.numpy_dtype.type, np.float32) else False
    constant_suffix = "f" if is_float else ""

    required_headers = get_headers(stream_collide_ast)

    if refinement_scaling:
        refinement_scaling_info = [(e0, e1, expression_to_code(e2, '', dtype=default_dtype)) for e0, e1, e2 in
                                   refinement_scaling.scaling_info]
        # append '_' to entries since they are used as members
        for i in range(len(refinement_scaling_info)):
            updated_entry = (refinement_scaling_info[i][0],
                             refinement_scaling_info[i][1].replace(refinement_scaling_info[i][1],
                                                                   refinement_scaling_info[i][1] + '_'),
                             refinement_scaling_info[i][2].replace(refinement_scaling_info[i][1],
                                                                   refinement_scaling_info[i][1] + '_'),
                             )
            refinement_scaling_info[i] = updated_entry
    else:
        refinement_scaling_info = None

    jinja_context = {
        'class_name': class_name,
        'stencil_name': stencil_name,
        'communication_stencil_name': communication_stencil_name,
        'D': lb_method.stencil.D,
        'Q': lb_method.stencil.Q,
        'compressible': lb_method.conserved_quantity_computation.compressible,
        'zero_centered': lb_method.conserved_quantity_computation.zero_centered_pdfs,
        'equilibrium_deviation_only': lb_method.equilibrium_distribution.deviation_only,
        'weights': ",".join(str(w.evalf()) + constant_suffix for w in lb_method.weights),
        'inverse_weights': ",".join(str((1 / w).evalf()) + constant_suffix for w in lb_method.weights),
        'dtype': "float" if is_float else "double",

        'equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, equilibrium),
        'symmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, symmetric_equilibrium),
        'asymmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, asymmetric_equilibrium),
        'equilibrium': [cpp_printer.doprint(e.rhs) for e in equilibrium],

        'macroscopic_velocity_shift': macroscopic_velocity_shift,
        'density_getters': equations_to_code(cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=default_dtype),
        'momentum_density_getter': equations_to_code(momentum_density_getter, variables_without_prefix=pdfs_sym,
                                                     dtype=default_dtype),
        'density_velocity_setter_macroscopic_values': density_velocity_setter_macroscopic_values,

        'refinement_scaling_info': refinement_scaling_info,

        'stream_collide_kernel': KernelInfo(stream_collide_ast, ['pdfs_tmp'], [('pdfs', 'pdfs_tmp')], []),
        'collide_kernel': KernelInfo(collide_ast, [], [], []),
        'stream_kernel': KernelInfo(stream_ast, ['pdfs_tmp'], [('pdfs', 'pdfs_tmp')], []),
        'target': 'cpu',
        'namespace': 'lbm',
        'headers': required_headers,
        'need_block_offsets': [
            'block_offset_{}'.format(i) in [param.symbol.name for param in stream_collide_ast.get_parameters()] for i in
            range(3)],
    }

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('LatticeModel.tmpl.h').render(**jinja_context)
    source = env.get_template('LatticeModel.tmpl.cpp').render(**jinja_context)

    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.cpp", source)


def generate_lattice_model(generation_context, class_name, collision_rule, field_layout='fzyx', refinement_scaling=None,
                           target=Target.CPU, data_type=None, cpu_openmp=None, cpu_vectorize_info=None,
                           **create_kernel_params):

    config = config_from_context(generation_context, target=target, data_type=data_type,
                                 cpu_openmp=cpu_openmp, cpu_vectorize_info=cpu_vectorize_info, **create_kernel_params)

    # usually a numpy layout is chosen by default i.e. xyzf - which is bad for waLBerla where at least the spatial
    # coordinates should be ordered in reverse direction i.e. zyx
    lb_method = collision_rule.method

    q = lb_method.stencil.Q
    dim = lb_method.stencil.D

    if config.target == Target.GPU:
        raise ValueError("Lattice Models can only be generated for CPUs. To generate LBM on GPUs use sweeps directly")

    if field_layout == 'fzyx':
        config.cpu_vectorize_info['assume_inner_stride_one'] = True
    elif field_layout == 'zyxf':
        config.cpu_vectorize_info['assume_inner_stride_one'] = False

    src_field = ps.Field.create_generic('pdfs', dim, config.data_type['pdfs'].numpy_dtype,
                                        index_dimensions=1, layout=field_layout, index_shape=(q,))
    dst_field = ps.Field.create_generic('pdfs_tmp', dim, config.data_type['pdfs_tmp'].numpy_dtype,
                                        index_dimensions=1, layout=field_layout,
                                        index_shape=(q,))

    stream_collide_update_rule = create_lbm_kernel(collision_rule, src_field, dst_field, StreamPullTwoFieldsAccessor())
    stream_collide_ast = create_kernel(stream_collide_update_rule, config=config)
    stream_collide_ast.function_name = 'kernel_streamCollide'
    stream_collide_ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']

    collide_update_rule = create_lbm_kernel(collision_rule, src_field, dst_field, CollideOnlyInplaceAccessor())
    collide_ast = create_kernel(collide_update_rule, config=config)
    collide_ast.function_name = 'kernel_collide'
    collide_ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']

    stream_update_rule = create_stream_only_kernel(lb_method.stencil, src_field, dst_field,
                                                   accessor=StreamPullTwoFieldsAccessor())
    stream_ast = create_kernel(stream_update_rule, config=config)
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = config.cpu_vectorize_info['assume_inner_stride_one']
    __lattice_model(generation_context, class_name, config, lb_method, stream_collide_ast, collide_ast, stream_ast,
                    refinement_scaling)


class RefinementScaling:
    level_scale_factor = sp.Symbol("level_scale_factor")

    def __init__(self):
        self.scaling_info = []

    def add_standard_relaxation_rate_scaling(self, viscosity_relaxation_rate):
        self.add_scaling(viscosity_relaxation_rate, relaxation_rate_scaling)

    def add_force_scaling(self, force_parameter):
        self.add_scaling(force_parameter, lambda param, factor: param * (1 / factor))

    def add_scaling(self, parameter, scaling_rule):
        """
        Adds a scaling rule, how parameters on refined blocks are modified

        :param parameter: parameter to modify: may either be a Field, Field.Access or a Symbol
        :param scaling_rule: function taking the parameter to be scaled as symbol and the scaling factor i.e.
                            how much finer the current block is compared to coarsest resolution
        """
        if isinstance(parameter, Field):
            field = parameter
            name = field.name
            if field.index_dimensions > 0:
                scaling_type = 'field_with_f'
                field_access = field(sp.Symbol("f"))
            else:
                scaling_type = 'field_xyz'
                field_access = field.center
            expr = scaling_rule(field_access, self.level_scale_factor)
            self.scaling_info.append((scaling_type, name, expr))
        elif isinstance(parameter, Field.Access):
            field_access = parameter
            expr = scaling_rule(field_access, self.level_scale_factor)
            name = field_access.field.name
            self.scaling_info.append(('field_xyz', name, expr))
        elif isinstance(parameter, sp.Symbol):
            expr = scaling_rule(parameter, self.level_scale_factor)
            self.scaling_info.append(('normal', parameter.name, expr))
        elif isinstance(parameter, list) or isinstance(parameter, tuple):
            for p in parameter:
                self.add_scaling(p, scaling_rule)
        else:
            raise ValueError("Invalid value for viscosity_relaxation_rate")


# ------------------------------------------ Internal ------------------------------------------------------------------


def stencil_switch_statement(stencil, values):
    template = Template("""
    using namespace stencil;
    switch( direction ) {
        {% for direction_name, value in dir_to_value_dict.items() -%}
            case {{direction_name}}: return {{value}};
        {% endfor -%}
        default:
            WALBERLA_ABORT("Invalid Direction")
    }
    """)

    dir_to_value_dict = {offset_to_direction_string(d): cpp_printer.doprint(v.rhs) for d, v in zip(stencil, values)}
    return template.render(dir_to_value_dict=dir_to_value_dict, undefined=StrictUndefined)


def field_and_symbol_substitute(expr, variable_prefix="lm.", variables_without_prefix=None):
    if variables_without_prefix is None:
        variables_without_prefix = []
    variables_without_prefix = [v.name if isinstance(v, sp.Symbol) else v for v in variables_without_prefix]
    substitutions = {}
    # check for member access
    if variable_prefix.endswith('.'):
        postfix = '_'
    else:
        postfix = ''
    for sym in expr.atoms(sp.Symbol):
        if isinstance(sym, Field.Access):
            fa = sym
            prefix = "" if fa.field.name in variables_without_prefix else variable_prefix
            if prefix.endswith('.'):
                postfix2 = '_'
            else:
                postfix2 = ''
            if fa.field.index_dimensions == 0:
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z)")
            else:
                assert fa.field.index_dimensions == 1, "walberla supports only 0 or 1 index dimensions"
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z,{fa.index[0]})")
        else:
            if sym.name not in variables_without_prefix:
                substitutions[sym] = sp.Symbol(variable_prefix + sym.name + postfix)
    return expr.subs(substitutions)


def expression_to_code(expr, variable_prefix="lm.", variables_without_prefix=None, dtype=None):
    """
    Takes a sympy expression and creates a C code string from it. Replaces field accesses by
    walberla field accesses i.e. field_W^1 -> field->get(-1, 0, 0, 1)
    :param dtype: default data type used for numbers in the code
    :param expr: sympy expression
    :param variable_prefix: all variables (and field) are prefixed with this string
                           this is used for member variables mostly
    :param variables_without_prefix: this variables are not prefixed
    :return: code string
    """
    if dtype is None:
        dtype = BasicType("float64")

    if variables_without_prefix is None:
        variables_without_prefix = []
    return cpp_printer.doprint(
        type_expr(field_and_symbol_substitute(expr, variable_prefix, variables_without_prefix), dtype=dtype))


def type_expr(eq, dtype):
    # manually cast 0.5 to dtype since this is somehow not done automatically
    eq = eq.subs(sp.Rational(1, 2), CastFunc(sp.Rational(1, 2), dtype))
    return eq.subs({s: TypedSymbol(s.name, dtype) for s in eq.atoms(sp.Symbol)})


def equations_to_code(equations, variable_prefix="lm.", variables_without_prefix=None, dtype=None):
    if dtype is None:
        dtype = BasicType("float64")

    if variables_without_prefix is None:
        variables_without_prefix = []
    if isinstance(equations, AssignmentCollection):
        equations = equations.all_assignments

    variables_without_prefix = list(variables_without_prefix)
    c_backend = CBackend()
    result = []
    left_hand_side_names = [e.lhs.name for e in equations]
    for eq in equations:
        assignment = SympyAssignment(type_expr(eq.lhs, dtype=dtype),
                                     type_expr(field_and_symbol_substitute(eq.rhs, variable_prefix,
                                                                           variables_without_prefix
                                                                           + left_hand_side_names),
                                               dtype=dtype))
        result.append(c_backend(assignment))
    return "\n".join(result)
