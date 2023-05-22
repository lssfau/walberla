from functools import reduce

from pystencils import Target

from pystencils.backends.cbackend import get_headers
from pystencils.backends.cuda_backend import CudaSympyPrinter
from pystencils.typing.typed_sympy import SHAPE_DTYPE
from pystencils.typing import TypedSymbol

from pystencils_walberla.utility import merge_sorted_lists


# TODO KernelInfo and KernelFamily should have same interface
class KernelInfo:
    def __init__(self, ast, temporary_fields=(), field_swaps=(), varying_parameters=()):
        self.ast = ast
        self.temporary_fields = tuple(temporary_fields)
        self.field_swaps = tuple(field_swaps)
        self.varying_parameters = tuple(varying_parameters)
        self.parameters = ast.get_parameters()  # cache parameters here

    @property
    def fields_accessed(self):
        return self.ast.fields_accessed

    def get_ast_attr(self, name):
        """Returns the value of an attribute of the AST managed by this KernelInfo.
        For compatibility with KernelFamily."""
        return self.ast.__getattribute__(name)

    def get_headers(self):
        all_headers = [list(get_headers(self.ast))]
        return reduce(merge_sorted_lists, all_headers)

    def generate_kernel_invocation_code(self, **kwargs):
        ast = self.ast
        ast_params = self.parameters
        is_cpu = self.ast.target == Target.CPU
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
