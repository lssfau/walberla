from os import path
from pystencils.data_types import get_base_type
from pystencils_walberla.cmake_integration import CodeGenerationContext

HEADER_EXTENSIONS = {'.h', '.hpp'}


def generate_info_header(ctx: CodeGenerationContext,
                         filename: str,
                         stencil_typedefs: dict = None,
                         field_typedefs: dict = None,
                         additional_headers: set = None,
                         headers_to_ignore: set = None,
                         additional_typedefs: dict = None,
                         additional_code: str = ""):
    """Generates an info header, consolidating required information about the generated code.
    The info header #includes all generated header files, and is thus the only header the
    application needs to #include. It can also contain aliases for waLBerla stencil types and
    instances of the GhostLayerField template.

    Args:
        ctx: Code Generation Context
        filename: Name of the generated info header file
        stencil_typedefs: dict mapping type names to stencil names or tuples
        field_typedefs: dict mapping type names to pystencils `Field` instances
        additional_headers: additional header files to be included
        headers_to_ignore: headers which should not be included
        additional_typedefs: dict mapping aliases to types
        additional_code: additional code which gets appended on the file
    """
    stencil_typedefs = stencil_typedefs if stencil_typedefs is not None else dict()
    field_typedefs = field_typedefs if field_typedefs is not None else dict()
    additional_typedefs = additional_typedefs if additional_typedefs is not None else dict()

    additional_headers = additional_headers if additional_headers is not None else set()
    headers_to_ignore = headers_to_ignore if headers_to_ignore is not None else set()

    headers_in_ctx = set(_filter_headers(ctx.files_written))

    stencil_headers, stencil_typedefs = _stencil_inclusion_code(stencil_typedefs)
    field_headers, field_typedefs = _field_inclusion_code(field_typedefs)

    headers_to_include = stencil_headers | field_headers | headers_in_ctx | additional_headers
    headers_to_include = sorted(headers_to_include - headers_to_ignore)
    headers_to_include = [f'#include "{header}"' for header in headers_to_include]

    typedefs = {**stencil_typedefs, **field_typedefs, **additional_typedefs}
    typedefs = [f"using {alias} = {typename};" for alias, typename in typedefs.items()]

    lines = '\n'.join(headers_to_include + [''] + typedefs) + '\n'

    if path.splitext(filename)[1] not in HEADER_EXTENSIONS:
        filename += '.h'

    ctx.write_file(filename, lines + additional_code)


#   ------------------------------------- INTERNAL -------------------------------------------------------------


def _filter_headers(filepaths):
    for p in filepaths:
        if path.splitext(p)[1] in HEADER_EXTENSIONS:
            yield path.split(p)[1]


def _stencil_inclusion_code(stencil_typedefs):
    headers = set()
    typedefs = dict()
    for typename, stencil in stencil_typedefs.items():
        if isinstance(stencil, tuple):
            dim = len(stencil[0])
            q = len(stencil)
            stencil = f"D{dim}Q{q}"
        elif not isinstance(stencil, str):
            raise ValueError(f'Invalid stencil: Do not know what to make of {stencil}')

        headers.add(f"stencil/{stencil}.h")
        typedefs[typename] = f"walberla::stencil::{stencil}"

    return headers, typedefs


def _field_inclusion_code(field_typedefs):
    headers = set()
    typedefs = dict()
    for typename, field in field_typedefs.items():
        f_size = field.values_per_cell()
        dtype = get_base_type(field.dtype)
        headers.add("field/GhostLayerField.h")
        typedefs[typename] = f"walberla::field::GhostLayerField<{dtype}, {f_size}>"

    return headers, typedefs
