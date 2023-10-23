# import warnings

from dataclasses import replace
from jinja2 import Environment, PackageLoader, StrictUndefined
import numpy as np

from pystencils import fields, Target

from lbmpy import LBMConfig, LBMOptimisation
from lbmpy.advanced_streaming import is_inplace, get_accessor, Timestep
from lbmpy.methods import AbstractLbMethod

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.utility import config_from_context
from lbmpy_walberla.packing_kernels import PackingKernelsCodegen
from lbmpy_walberla.utility import create_pdf_field


def generate_lbm_storage_specification(generation_context: CodeGenerationContext, class_name: str,
                                       method: AbstractLbMethod,
                                       lbm_config: LBMConfig, lbm_optimisation: LBMOptimisation,
                                       nonuniform: bool = False,
                                       target: Target = Target.CPU,
                                       data_type=None, cpu_openmp: bool = False,
                                       **create_kernel_params):
    namespace = "lbm"
    stencil = method.stencil
    stencil_name = stencil.name
    streaming_pattern = lbm_config.streaming_pattern

    config = config_from_context(generation_context, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 **create_kernel_params)

    # Packing kernels should never be vectorised
    config = replace(config, cpu_vectorize_info=None)

    default_dtype = config.data_type.default_factory()
    if issubclass(default_dtype.numpy_dtype.type, np.float64):
        data_type_string = "double"
    elif issubclass(default_dtype.numpy_dtype.type, np.float32):
        data_type_string = "float"
    elif issubclass(default_dtype.numpy_dtype.type, np.float16):
        data_type_string = "half"
    else:
        raise ValueError(f"default datatype {default_dtype.numpy_dtype.type} is not supported. "
                         f"Supported are only np.float64, np.float32 and np.float16")

    symbolic_field = lbm_optimisation.symbolic_field
    if not symbolic_field:
        symbolic_field = create_pdf_field(config=config, name="pdfs_src", stencil=stencil,
                                          field_layout=lbm_optimisation.field_layout)

    if is_inplace(streaming_pattern):
        symbolic_temporary_field = create_pdf_field(config=config, name="pdfs_dst", stencil=stencil,
                                                    field_layout=lbm_optimisation.field_layout)
    else:
        symbolic_temporary_field = lbm_optimisation.symbolic_temporary_field
        if not symbolic_temporary_field:
            symbolic_temporary_field = create_pdf_field(config=config, name="pdfs_dst", stencil=stencil,
                                                        field_layout=lbm_optimisation.field_layout)

    cg = PackingKernelsCodegen(stencil, streaming_pattern, class_name, config,
                               src_field=symbolic_field, dst_field=symbolic_temporary_field)
    kernels = cg.create_uniform_kernel_families()

    if nonuniform:
        kernels = cg.create_nonuniform_kernel_families(kernels_dict=kernels)

    # Pure storage specification
    if not stencil_name:
        raise ValueError("lb_method uses a stencil that is not supported in waLBerla")

    communication_stencil_name = stencil_name if stencil_name != "D3Q15" else "D3Q27"

    cqc = method.conserved_quantity_computation
    equilibrium = method.equilibrium_distribution

    f = fields(f"f({stencil.Q}): double[{stencil.D}D]", layout='fzyx')
    even_accessor = get_accessor(streaming_pattern, Timestep.EVEN)
    odd_accessor = get_accessor(streaming_pattern, Timestep.ODD)

    even_read = even_accessor.read(f, stencil)
    even_write = even_accessor.write(f, stencil)

    odd_read = odd_accessor.read(f, stencil)
    odd_write = odd_accessor.write(f, stencil)

    jinja_context = {
        'class_name': class_name,
        'namespace': namespace,
        'stencil_name': stencil_name,
        'communication_stencil_name': communication_stencil_name,
        'stencil_size': stencil.Q,
        'dimension': stencil.D,
        'compressible': cqc.compressible,
        'equilibrium_accuracy_order': equilibrium.order,
        'equilibrium_deviation_only': equilibrium.deviation_only,
        'inplace': is_inplace(streaming_pattern),
        'zero_centered': cqc.zero_centered_pdfs,

        'weights': ", ".join(f"{data_type_string}({str(w.evalf())})" for w in method.weights),
        'inverse_weights': ", ".join(f"{data_type_string}({str((1 / w).evalf())})" for w in method.weights),
        'even_read': _get_access_list(even_read, stencil.D),
        'even_write': _get_access_list(even_write, stencil.D),
        'odd_read': _get_access_list(odd_read, stencil.D),
        'odd_write': _get_access_list(odd_write, stencil.D),

        'nonuniform': nonuniform,
        'target': target.name.lower(),
        'dtype': data_type_string,
        'is_gpu': target == Target.GPU,
        'kernels': kernels,
        'direction_sizes': cg.get_direction_sizes(),
        'src_field': cg.src_field,
        'dst_field': cg.dst_field

    }
    if nonuniform:
        jinja_context['mask_field'] = cg.mask_field

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('LbmStorageSpecification.tmpl.h').render(**jinja_context)
    source = env.get_template('LbmStorageSpecification.tmpl.cpp').render(**jinja_context)

    source_extension = "cu" if target == Target.GPU and generation_context.cuda else "cpp"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)


def _get_access_list(access_list, dim):
    result = []
    for i in range(dim):
        result.append(", ".join([str(int(field_access.offsets[i])) for field_access in access_list]))

    if dim == 2:
        result.append(", ".join(["0"] * len(access_list)))

    result.append(", ".join([str(int(field_access.index[0])) for field_access in access_list]))

    return result
