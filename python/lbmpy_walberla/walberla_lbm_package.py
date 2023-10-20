from typing import Callable, List, Dict

from pystencils import Target, Field

from lbmpy.creationfunctions import LbmCollisionRule, LBMConfig, LBMOptimisation
from lbmpy.relaxationrates import get_shear_relaxation_rate

from pystencils_walberla.cmake_integration import CodeGenerationContext

from lbmpy_walberla.boundary_collection import generate_boundary_collection
from lbmpy_walberla.storage_specification import generate_lbm_storage_specification
from lbmpy_walberla.sweep_collection import generate_lbm_sweep_collection, RefinementScaling


def generate_lbm_package(ctx: CodeGenerationContext, name: str,
                         collision_rule: LbmCollisionRule,
                         lbm_config: LBMConfig, lbm_optimisation: LBMOptimisation,
                         nonuniform: bool = False, boundaries: List[Callable] = None,
                         macroscopic_fields: Dict[str, Field] = None,
                         target: Target = Target.CPU,
                         data_type=None, pdfs_data_type=None,
                         cpu_openmp=None, cpu_vectorize_info=None,
                         max_threads=None,
                         **kernel_parameters):

    if macroscopic_fields is None:
        macroscopic_fields = {}

    method = collision_rule.method

    storage_spec_name = f'{name}StorageSpecification'
    generate_lbm_storage_specification(ctx, storage_spec_name, method, lbm_config, lbm_optimisation,
                                       nonuniform=nonuniform, target=target,
                                       data_type=pdfs_data_type,
                                       cpu_openmp=cpu_openmp)

    if nonuniform:
        omega = get_shear_relaxation_rate(method)
        refinement_scaling = RefinementScaling()
        refinement_scaling.add_standard_relaxation_rate_scaling(omega)
    else:
        refinement_scaling = None

    generate_lbm_sweep_collection(ctx, f'{name}SweepCollection', collision_rule,
                                  lbm_config=lbm_config, lbm_optimisation=lbm_optimisation,
                                  refinement_scaling=refinement_scaling,
                                  macroscopic_fields=macroscopic_fields,
                                  target=target, data_type=data_type,
                                  cpu_openmp=cpu_openmp, cpu_vectorize_info=cpu_vectorize_info,
                                  max_threads=max_threads,
                                  **kernel_parameters)

    spatial_shape = None
    if lbm_optimisation.symbolic_field and lbm_optimisation.symbolic_field.has_fixed_shape:
        spatial_shape = lbm_optimisation.symbolic_field.spatial_shape + (lbm_config.stencil.Q, )

    generate_boundary_collection(ctx, f'{name}BoundaryCollection', boundary_generators=boundaries,
                                 lb_method=method, field_name='pdfs', spatial_shape=spatial_shape,
                                 streaming_pattern=lbm_config.streaming_pattern,
                                 target=target, layout=lbm_optimisation.field_layout)
