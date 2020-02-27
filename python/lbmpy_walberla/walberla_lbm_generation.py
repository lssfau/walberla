import warnings

import numpy as np
import sympy as sp
from jinja2 import Environment, PackageLoader, StrictUndefined, Template
from sympy.tensor import IndexedBase

import pystencils as ps
from lbmpy.creationfunctions import create_lb_update_rule, update_with_default_parameters
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor, StreamPullTwoFieldsAccessor
from lbmpy.relaxationrates import relaxation_rate_scaling
from lbmpy.stencils import get_stencil
from lbmpy.updatekernels import create_lbm_kernel, create_stream_pull_only_kernel
from pystencils import AssignmentCollection, create_kernel
from pystencils.astnodes import SympyAssignment
from pystencils.backends.cbackend import CBackend, CustomSympyPrinter, get_headers
from pystencils.data_types import TypedSymbol, type_all_numbers
from pystencils.field import Field
from pystencils.stencil import have_same_entries, offset_to_direction_string
from pystencils.sympyextensions import get_symmetric_part
from pystencils.transformations import add_types
from pystencils_walberla.codegen import KernelInfo, default_create_kernel_parameters
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env

cpp_printer = CustomSympyPrinter()
REFINEMENT_SCALE_FACTOR = sp.Symbol("level_scale_factor")


def __lattice_model(generation_context, class_name, lb_method, stream_collide_ast, collide_ast, stream_ast,
                    refinement_scaling):
    stencil_name = get_stencil_name(lb_method.stencil)
    if not stencil_name:
        raise ValueError("lb_method uses a stencil that is not supported in waLBerla")
    is_float = not generation_context.double_accuracy
    dtype = np.float32 if is_float else np.float64

    vel_symbols = lb_method.conserved_quantity_computation.first_order_moment_symbols
    rho_sym = sp.Symbol("rho")
    pdfs_sym = sp.symbols("f_:%d" % (len(lb_method.stencil),))
    vel_arr_symbols = [IndexedBase(sp.Symbol('u'), shape=(1,))[i] for i in range(len(vel_symbols))]
    momentum_density_symbols = [sp.Symbol("md_%d" % (i,)) for i in range(len(vel_symbols))]

    equilibrium = lb_method.get_equilibrium()
    equilibrium = equilibrium.new_with_substitutions({a: b for a, b in zip(vel_symbols, vel_arr_symbols)})
    _, _, equilibrium = add_types(equilibrium.main_assignments, "float32" if is_float else "float64", False)
    equilibrium = sp.Matrix([e.rhs for e in equilibrium])

    symmetric_equilibrium = get_symmetric_part(equilibrium, vel_arr_symbols)
    asymmetric_equilibrium = sp.expand(equilibrium - symmetric_equilibrium)

    force_model = lb_method.force_model
    macroscopic_velocity_shift = None
    if force_model:
        if hasattr(force_model, 'macroscopic_velocity_shift'):
            macroscopic_velocity_shift = [expression_to_code(e, "lm.", ["rho"],dtype=dtype)
                                          for e in force_model.macroscopic_velocity_shift(rho_sym)]

    cqc = lb_method.conserved_quantity_computation

    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(sp.Symbol("rho_in"), vel_arr_symbols)
    density_velocity_setter_macroscopic_values = equations_to_code(eq_input_from_input_eqs, dtype=dtype,
                                                                   variables_without_prefix=['rho_in', 'u'])
    momentum_density_getter = cqc.output_equations_from_pdfs(pdfs_sym, {'density': rho_sym,
                                                                        'momentum_density': momentum_density_symbols})
    constant_suffix = "f" if is_float else ""

    required_headers = get_headers(stream_collide_ast)

    if refinement_scaling:
        refinement_scaling_info = [ (e0,e1,expression_to_code(e2, '', dtype=dtype)) for e0,e1,e2 in refinement_scaling.scaling_info ]
    else:
        refinement_scaling_info = None

    jinja_context = {
        'class_name': class_name,
        'stencil_name': stencil_name,
        'D': lb_method.dim,
        'Q': len(lb_method.stencil),
        'compressible': lb_method.conserved_quantity_computation.compressible,
        'weights': ",".join(str(w.evalf()) + constant_suffix for w in lb_method.weights),
        'inverse_weights': ",".join(str((1/w).evalf()) + constant_suffix for w in lb_method.weights),

        'equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, equilibrium),
        'symmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, symmetric_equilibrium),
        'asymmetric_equilibrium_from_direction': stencil_switch_statement(lb_method.stencil, asymmetric_equilibrium),
        'equilibrium': [cpp_printer.doprint(e) for e in equilibrium],

        'macroscopic_velocity_shift': macroscopic_velocity_shift,
        'density_getters': equations_to_code(cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=dtype),
        'momentum_density_getter': equations_to_code(momentum_density_getter, variables_without_prefix=pdfs_sym,
                                                     dtype=dtype),
        'density_velocity_setter_macroscopic_values': density_velocity_setter_macroscopic_values,

        'refinement_scaling_info': refinement_scaling_info,

        'stream_collide_kernel': KernelInfo(stream_collide_ast, ['pdfs_tmp'], [('pdfs', 'pdfs_tmp')], []),
        'collide_kernel': KernelInfo(collide_ast, [], [], []),
        'stream_kernel': KernelInfo(stream_ast, ['pdfs_tmp'], [('pdfs', 'pdfs_tmp')], []),
        'target': 'cpu',
        'namespace': 'lbm',
        'headers': required_headers,
        'need_block_offsets': ['block_offset_{}'.format(i) in [param.symbol.name for param in stream_collide_ast.get_parameters()] for i in range(3)],
    }

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('LatticeModel.tmpl.h').render(**jinja_context)
    source = env.get_template('LatticeModel.tmpl.cpp').render(**jinja_context)

    generation_context.write_file("{}.h".format(class_name), header)
    generation_context.write_file("{}.cpp".format(class_name), source)


def generate_lattice_model(generation_context, class_name, collision_rule, refinement_scaling=None,
                           **create_kernel_params):

    # usually a numpy layout is chosen by default i.e. xyzf - which is bad for waLBerla where at least the spatial
    # coordinates should be ordered in reverse direction i.e. zyx
    is_float = not generation_context.double_accuracy
    dtype = np.float32 if is_float else np.float64
    lb_method = collision_rule.method

    q = len(lb_method.stencil)
    dim = lb_method.dim

    create_kernel_params = default_create_kernel_parameters(generation_context, create_kernel_params)
    if create_kernel_params['target'] == 'gpu':
        raise ValueError("Lattice Models can only be generated for CPUs. To generate LBM on GPUs use sweeps directly")

    src_field = ps.Field.create_generic('pdfs', dim, dtype, index_dimensions=1, layout='fzyx', index_shape=(q,))
    dst_field = ps.Field.create_generic('pdfs_tmp', dim, dtype, index_dimensions=1, layout='fzyx', index_shape=(q,))

    stream_collide_update_rule = create_lbm_kernel(collision_rule, src_field, dst_field, StreamPullTwoFieldsAccessor())
    stream_collide_ast = create_kernel(stream_collide_update_rule, **create_kernel_params)
    stream_collide_ast.function_name = 'kernel_streamCollide'

    collide_update_rule = create_lbm_kernel(collision_rule, src_field, dst_field, CollideOnlyInplaceAccessor())
    collide_ast = create_kernel(collide_update_rule, **create_kernel_params)
    collide_ast.function_name = 'kernel_collide'

    stream_update_rule = create_stream_pull_only_kernel(lb_method.stencil, None, 'pdfs', 'pdfs_tmp', 'fzyx', dtype)
    stream_ast = create_kernel(stream_update_rule, **create_kernel_params)
    stream_ast.function_name = 'kernel_stream'
    __lattice_model(generation_context, class_name, lb_method, stream_collide_ast, collide_ast, stream_ast,
                    refinement_scaling)




class RefinementScaling:
    level_scale_factor = sp.Symbol("level_scale_factor")

    def __init__(self):
        self.scaling_info = []

    def add_standard_relaxation_rate_scaling(self, viscosity_relaxation_rate):
        self.add_scaling(viscosity_relaxation_rate, relaxation_rate_scaling)

    def add_force_scaling(self, force_parameter):
        self.add_scaling(force_parameter, lambda param, factor: param * factor)

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
            WALBERLA_ABORT("Invalid Direction");
    }
    """)

    dir_to_value_dict = {offset_to_direction_string(d): cpp_printer.doprint(v) for d, v in zip(stencil, values)}
    return template.render(dir_to_value_dict=dir_to_value_dict, undefined=StrictUndefined)


def field_and_symbol_substitute(expr, variable_prefix="lm.", variables_without_prefix=[]):
    variables_without_prefix = [v.name if isinstance(v, sp.Symbol) else v for v in variables_without_prefix]
    substitutions = {}
    for sym in expr.atoms(sp.Symbol):
        if isinstance(sym, Field.Access):
            fa = sym
            prefix = "" if fa.field.name in variables_without_prefix else variable_prefix
            if fa.field.index_dimensions == 0:
                substitutions[fa] = sp.Symbol("%s%s->get(x,y,z)" % (prefix, fa.field.name))
            else:
                assert fa.field.index_dimensions == 1, "walberla supports only 0 or 1 index dimensions"
                substitutions[fa] = sp.Symbol("%s%s->get(x,y,z,%s)" % (prefix, fa.field.name, fa.index[0]))
        else:
            if sym.name not in variables_without_prefix:
                substitutions[sym] = sp.Symbol(variable_prefix + sym.name)
    return expr.subs(substitutions)


def expression_to_code(expr, variable_prefix="lm.", variables_without_prefix=[],dtype="double"):
    """
    Takes a sympy expression and creates a C code string from it. Replaces field accesses by
    walberla field accesses i.e. field_W^1 -> field->get(-1, 0, 0, 1)
    :param expr: sympy expression
    :param variable_prefix: all variables (and field) are prefixed with this string
                           this is used for member variables mostly
    :param variables_without_prefix: this variables are not prefixed
    :return: code string
    """
    return cpp_printer.doprint(type_expr(field_and_symbol_substitute(expr, variable_prefix, variables_without_prefix),dtype=dtype))

def type_expr(eq, dtype):
    eq=type_all_numbers(eq,dtype=dtype)
    return eq.subs({s: TypedSymbol(s.name, dtype) for s in eq.atoms(sp.Symbol)})

def equations_to_code(equations, variable_prefix="lm.", variables_without_prefix=[], dtype="double"):

    if isinstance(equations, AssignmentCollection):
        equations = equations.all_assignments

    variables_without_prefix = list(variables_without_prefix)
    c_backend = CBackend()
    result = []
    left_hand_side_names = [e.lhs.name for e in equations]
    for eq in equations:
        assignment = SympyAssignment(type_expr(eq.lhs,dtype=dtype),
                                     type_expr(field_and_symbol_substitute(eq.rhs, variable_prefix,
                                                                 variables_without_prefix + left_hand_side_names),dtype=dtype))
        result.append(c_backend(assignment))
    return "\n".join(result)


def get_stencil_name(stencil):
    for name in ('D2Q9', 'D3Q15', 'D3Q19', 'D3Q27'):
        if have_same_entries(stencil, get_stencil(name, 'walberla')):
            return name
