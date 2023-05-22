# import warnings

from dataclasses import replace
from jinja2 import Environment, PackageLoader, StrictUndefined
import numpy as np

from pystencils import Target

from lbmpy import LBMConfig
from lbmpy.advanced_streaming import is_inplace
from lbmpy.methods import AbstractLbMethod

from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.utility import config_from_context
from lbmpy_walberla.packing_kernels import PackingKernelsCodegen


def generate_lbm_storage_specification(generation_context, class_name: str,
                                       method: AbstractLbMethod, lbm_config: LBMConfig, nonuniform: bool = False,
                                       target: Target = Target.CPU, data_type=None, cpu_openmp: bool = False,
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
    is_float = True if issubclass(default_dtype.numpy_dtype.type, np.float32) else False

    cg = PackingKernelsCodegen(stencil, streaming_pattern, class_name, config)
    kernels = cg.create_uniform_kernel_families()

    if nonuniform:
        kernels = cg.create_nonuniform_kernel_families(kernels_dict=kernels)

    values_per_cell = len(stencil)
    dimension = len(stencil[0])

    # Pure storage specification
    if not stencil_name:
        raise ValueError("lb_method uses a stencil that is not supported in waLBerla")

    communication_stencil_name = stencil_name if stencil_name != "D3Q15" else "D3Q27"

    cqc = method.conserved_quantity_computation
    equilibrium = method.equilibrium_distribution

    jinja_context = {
        'class_name': class_name,
        'namespace': namespace,
        'stencil_name': stencil_name,
        'communication_stencil_name': communication_stencil_name,
        'compressible': cqc.compressible,
        'equilibriumAccuracyOrder': equilibrium.order,
        'inplace': is_inplace(streaming_pattern),
        'zero_centered': cqc.zero_centered_pdfs,
        'eq_deviation_only': equilibrium.deviation_only,

        'nonuniform': nonuniform,
        'target': target.name.lower(),
        'dtype': "float" if is_float else "double",
        'is_gpu': target == Target.GPU,
        'kernels': kernels,
        'direction_sizes': cg.get_direction_sizes(),
        'stencil_size': stencil.Q,
        'dimension': stencil.D,
        'src_field': cg.src_field,
        'dst_field': cg.dst_field

    }
    if nonuniform:
        jinja_context['mask_field'] = cg.mask_field

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)

    header = env.get_template('LbmStorageSpecification.tmpl.h').render(**jinja_context)
    source = env.get_template('LbmStorageSpecification.tmpl.cpp').render(**jinja_context)

    source_extension = "cpp" if target == Target.CPU else "cu"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)
