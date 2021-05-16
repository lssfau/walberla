from typing import Sequence, Set
from collections import OrderedDict
from functools import reduce
from jinja2.filters import do_indent
from pystencils import TypedSymbol
from pystencils.backends.cbackend import get_headers
from pystencils.backends.cuda_backend import CudaSympyPrinter
from pystencils.kernelparameters import SHAPE_DTYPE


"""

This module contains several classes and methods supporting the code generation of sweeps
containing multiple kernels.

A sweep class may be backed by multiple kernel implementations with (nearly) the same interface.
Those might be slightly altered versions of a single kernel, and the one to be executed may
be selected at run time according to certain criteria. An example are the even/odd alternating
kernels required for lattice Boltzmann inplace streaming. Executing multiple kernels in sequence
within a single sweep could also be a use case, which is however not yet implemented.

## The Kernel Selection Tree

The selection of the correct kernel is modelled by a tree structure, spanned by instances of
subclasses of `AbstractKernelSelectionNode`. The code generator traverses this tree and generates
code for invoking kernels according to its structure and nodes.
Currently, two types of nodes exist:

- Condition Nodes: Subclasses of `AbstractConditionNode` manifest as if/else statements
    in the generated code. They model decision points and have two subtrees; one for the
    `true` case, and one for `false`. A basic implementation is `SimpleBooleanCondition`,
    which simply branches according to a boolean value. Each condition node requires a number
    of selection parameters. The total of selection parameters of the entire tree are collected
    during code generation, and must all be passed to the generated C++ functions.
- Kernel Call Nodes: Currently, `KernelCallNode` corresponds to the invocation of a single
    kernel. `KernelCallNode` acts as a wrapper for a single pystencils AST. When encountered
    in the tree during code generation, the function call for this AST is inserted.


## The Kernel Family

The `KernelFamily` class is a wrapper around a kernel selection tree. It was developed as a
generalization of the `KernelInfo` class. Its purpose is the collection and management of
information about the tree and its kernels which is required for code generation.
It also checks the tree's ASTs for consistency; for example by making sure that any fields
and symbols required by multiple ASTs have the same type, dimensions, et cetera.


## High-Level Interface

Due to the tree's selection arguments, which must be passed to the methods wrapping the
kernel calls, the generated class can not by itself be used as a sweep functor to be passed
to the waLBerla timeloop. Instead, for all sweep types, a `get[...]Sweep` member function is
generated. It takes any required selection arguments, and returns a lambda function which
can be passed directly to the timeloop.

Using the interface mapping system, the 'low-level' selection arguments of the kernel tree
can be hidden behind higher-level arguments. Such argument mappings are modelled by the
subclasses of `AbstractInterfaceArgumentMapping`. During code generation, they are organized
in an instance of `HighLevelInterfaceSpec`, which is used to generate the high-level interface.


"""

# ---------------------------------- Selection Tree --------------------------------------------------------------------


class AbstractKernelSelectionNode:

    def collect_kernel_calls(self):
        raise NotImplementedError()

    @property
    def selection_parameters(self) -> Set[TypedSymbol]:
        raise NotImplementedError()

    def collect_selection_parameters(self) -> Set[TypedSymbol]:
        return self.selection_parameters

    def get_selection_parameter_list(self) -> Sequence[TypedSymbol]:
        all_params = self.collect_selection_parameters()
        all_names = set(p.name for p in all_params)
        if len(all_names) < len(all_params):
            raise ValueError('There existed selection parameters of same name, but different type.')
        return sorted(all_params, key=lambda x: x.name)

    def get_code(self, **kwargs) -> str:
        raise NotImplementedError()


class AbstractConditionNode(AbstractKernelSelectionNode):
    def __init__(self, branch_true, branch_false):
        self.branch_true = branch_true
        self.branch_false = branch_false

    @property
    def condition_text(self) -> str:
        raise NotImplementedError()

    def collect_kernel_calls(self):
        return self.branch_true.collect_kernel_calls() | self.branch_false.collect_kernel_calls()

    def collect_selection_parameters(self) -> Set[TypedSymbol]:
        return self.selection_parameters.union(self.branch_true.collect_selection_parameters(),
                                               self.branch_false.collect_selection_parameters())

    def get_code(self, **kwargs):
        true_branch_code = self.branch_true.get_code(**kwargs)
        false_branch_code = self.branch_false.get_code(**kwargs)

        true_branch_code = do_indent(true_branch_code, width=4, first=True)
        false_branch_code = do_indent(false_branch_code, width=4, first=True)

        code = f"if({self.condition_text}) {{\n"
        code += true_branch_code
        code += "\n} else {\n"
        code += false_branch_code
        code += "\n}"
        return code


class KernelCallNode(AbstractKernelSelectionNode):
    def __init__(self, ast):
        self.ast = ast
        self.parameters = ast.get_parameters()  # cache parameters here

    @property
    def selection_parameters(self) -> Set[TypedSymbol]:
        return set()

    def collect_kernel_calls(self):
        return {self}

    def get_code(self, **kwargs):
        ast = self.ast
        ast_params = self.parameters
        is_cpu = self.ast.target == 'cpu'
        call_parameters = ", ".join([p.symbol.name for p in ast_params])

        if not is_cpu:
            stream = kwargs.get('stream', '0')
            spatial_shape_symbols = kwargs.get('spatial_shape_symbols', ())

            if not spatial_shape_symbols:
                spatial_shape_symbols = [p.symbol for p in ast_params if p.is_field_shape]
                spatial_shape_symbols.sort(key=lambda e: e.coordinate)
            else:
                spatial_shape_symbols = [TypedSymbol(s, SHAPE_DTYPE) for s in spatial_shape_symbols]

            assert spatial_shape_symbols, "No shape parameters in kernel function arguments.\n"\
                "Please only use kernels for generic field sizes!"

            indexing_dict = ast.indexing.call_parameters(spatial_shape_symbols)
            sp_printer_c = CudaSympyPrinter()
            kernel_call_lines = [
                "dim3 _block(int(%s), int(%s), int(%s));" % tuple(sp_printer_c.doprint(e)
                                                                  for e in indexing_dict['block']),
                "dim3 _grid(int(%s), int(%s), int(%s));" % tuple(sp_printer_c.doprint(e)
                                                                 for e in indexing_dict['grid']),
                "internal_%s::%s<<<_grid, _block, 0, %s>>>(%s);" % (ast.function_name, ast.function_name,
                                                                    stream, call_parameters),
            ]

            return "\n".join(kernel_call_lines)
        else:
            return f"internal_{ast.function_name}::{ast.function_name}({call_parameters});"


class SimpleBooleanCondition(AbstractConditionNode):
    def __init__(self,
                 parameter_name: str,
                 branch_true: AbstractKernelSelectionNode,
                 branch_false: AbstractKernelSelectionNode):
        self.parameter_symbol = TypedSymbol(parameter_name, bool)
        super(SimpleBooleanCondition, self).__init__(branch_true, branch_false)

    @property
    def selection_parameters(self) -> Set[TypedSymbol]:
        return {self.parameter_symbol}

    @property
    def condition_text(self) -> str:
        return self.parameter_symbol.name


# ---------------------------------- Kernel Family ---------------------------------------------------------------------


class KernelFamily:
    def __init__(self, kernel_selection_tree: AbstractKernelSelectionNode,
                 class_name: str,
                 temporary_fields=(), field_swaps=(), varying_parameters=()):
        self.kernel_selection_tree = kernel_selection_tree
        self.kernel_selection_parameters = kernel_selection_tree.get_selection_parameter_list()
        self.temporary_fields = tuple(temporary_fields)
        self.field_swaps = tuple(field_swaps)
        self.varying_parameters = tuple(varying_parameters)

        all_kernel_calls = self.kernel_selection_tree.collect_kernel_calls()
        all_param_lists = [k.parameters for k in all_kernel_calls]
        asts_list = [k.ast for k in all_kernel_calls]
        self.representative_ast = asts_list[0]

        #   Eliminate duplicates
        self.all_asts = set(asts_list)

        #   Check function names for uniqueness and reformat them
        #   using the class name
        function_names = [ast.function_name.lower() for ast in self.all_asts]
        unique_names = set(function_names)
        if len(unique_names) < len(function_names):
            raise ValueError('Function names of kernel family members must be unique!')

        prefix = class_name.lower()
        for ast in self.all_asts:
            ast.function_name = prefix + '_' + ast.function_name

        all_fields = [k.ast.fields_accessed for k in all_kernel_calls]

        #   Collect function parameters and accessed fields
        self.parameters = merge_lists_of_symbols(all_param_lists)
        self.fields_accessed = reduce(lambda x, y: x | y, all_fields)

        self._ast_attrs = dict()

    def get_ast_attr(self, name):
        """Returns the value of an attribute of the ASTs managed by this KernelFamily only
        if it is the same in all ASTs."""
        if name not in self._ast_attrs:
            attr = self.representative_ast.__getattribute__(name)
            for ast in self.all_asts:
                if ast.__getattribute__(name) != attr:
                    raise ValueError(f'Inconsistency in kernel family: Attribute {name} was different in {ast}!')
            self._ast_attrs[name] = attr
        return self._ast_attrs[name]

    def get_headers(self):
        all_headers = [get_headers(ast) for ast in self.all_asts]
        return reduce(merge_sorted_lists, all_headers)

    def generate_kernel_invocation_code(self, **kwargs):
        return self.kernel_selection_tree.get_code(**kwargs)


# --------------------------- High-Level Sweep Interface Specification ------------------------------------------------


class AbstractInterfaceArgumentMapping:
    def __init__(self, high_level_args: Sequence[TypedSymbol], low_level_arg: TypedSymbol):
        self.high_level_args = high_level_args
        self.low_level_arg = low_level_arg

    @property
    def mapping_code(self):
        raise NotImplementedError()

    @property
    def headers(self):
        return set()


class IdentityMapping(AbstractInterfaceArgumentMapping):

    def __init__(self, low_level_arg: TypedSymbol):
        self.high_level_args = (low_level_arg,)
        self.low_level_arg = low_level_arg

    @property
    def mapping_code(self):
        return self.low_level_arg.name


def _create_interface_mapping_dict(low_level_args: Sequence[TypedSymbol],
                                   mappings: Sequence[AbstractInterfaceArgumentMapping]):
    mapping_dict = OrderedDict()
    for m in mappings:
        larg = m.low_level_arg
        if larg not in low_level_args:
            raise ValueError(f'Low-level argument {larg} did not exist.')
        if larg.name in mapping_dict:
            raise ValueError(f'More than one mapping was given for low-level argument {larg}')
        mapping_dict[larg.name] = m

    for arg in low_level_args:
        mapping_dict.setdefault(arg.name, IdentityMapping(arg))

    return mapping_dict


class HighLevelInterfaceSpec:
    def __init__(self, low_level_args: Sequence[TypedSymbol],
                 mappings: Sequence[AbstractInterfaceArgumentMapping]):
        self.low_level_args = low_level_args
        mapping_dict = _create_interface_mapping_dict(low_level_args, mappings)
        self.mappings = mapping_dict.values()
        high_level_args = OrderedDict()
        self.mapping_codes = []
        self.headers = set()
        for larg in low_level_args:
            m = mapping_dict[larg.name]
            self.mapping_codes.append(m.mapping_code)
            self.headers |= m.headers
            for harg in m.high_level_args:
                if high_level_args.setdefault(harg.name, harg) != harg:
                    raise ValueError(f'High-Level Argument {harg} was given multiple times with different types.')

        self.high_level_args = list(high_level_args.values())


# ---------------------------------- Helpers --------------------------------------------------------------------------


def merge_sorted_lists(lx, ly, sort_key=lambda x: x, identity_check_key=None):
    if identity_check_key is None:
        identity_check_key = sort_key
    nx = len(lx)
    ny = len(ly)

    def recursive_merge(lx, ly, ix, iy):
        if ix == nx:
            return ly[iy:]
        if iy == ny:
            return lx[ix:]
        x = lx[ix]
        y = ly[iy]
        skx = sort_key(x)
        sky = sort_key(y)
        if skx == sky:
            if identity_check_key(x) == identity_check_key(y):
                return [x] + recursive_merge(lx, ly, ix + 1, iy + 1)
            else:
                raise ValueError(f'Elements <{x}> and <{y}> with equal sort key where not identical!')
        elif skx < sky:
            return [x] + recursive_merge(lx, ly, ix + 1, iy)
        else:
            return [y] + recursive_merge(lx, ly, ix, iy + 1)
    return recursive_merge(lx, ly, 0, 0)


def merge_lists_of_symbols(lists):
    def merger(lx, ly):
        return merge_sorted_lists(lx, ly, sort_key=lambda x: x.symbol.name, identity_check_key=lambda x: x.symbol)
    return reduce(merger, lists)
