from dataclasses import replace

import numpy as np

from pystencils_walberla.codegen import generate_selective_sweep, config_from_context
from pystencils_walberla.kernel_selection import (
    AbstractInterfaceArgumentMapping, AbstractConditionNode, KernelCallNode)
from pystencils import Target, TypedSymbol
from lbmpy.creationfunctions import create_lb_ast
from lbmpy.advanced_streaming import Timestep, is_inplace


class EvenIntegerCondition(AbstractConditionNode):
    def __init__(self, parameter_name: str,
                 branch_true, branch_false,
                 parameter_dtype=int):
        self.parameter_symbol = TypedSymbol(parameter_name, parameter_dtype)
        super(EvenIntegerCondition, self).__init__(branch_true, branch_false)

    @property
    def selection_parameters(self):
        return {self.parameter_symbol}

    @property
    def condition_text(self):
        return f"(({self.parameter_symbol.name} & 1) ^ 1)"


class OddIntegerCondition(AbstractConditionNode):
    def __init__(self, parameter_name: str,
                 branch_true, branch_false,
                 parameter_dtype=int):
        self.parameter_symbol = TypedSymbol(parameter_name, parameter_dtype)
        super(OddIntegerCondition, self).__init__(branch_true, branch_false)

    @property
    def selection_parameters(self):
        return {self.parameter_symbol}

    @property
    def condition_text(self):
        return f"({self.parameter_symbol.name} & 1)"


class TimestepTrackerMapping(AbstractInterfaceArgumentMapping):

    def __init__(self, low_level_arg: TypedSymbol, tracker_identifier='tracker'):
        self.tracker_symbol = TypedSymbol(tracker_identifier, 'std::shared_ptr<lbm::TimestepTracker> &')
        super(TimestepTrackerMapping, self).__init__(high_level_args=(self.tracker_symbol,),
                                                     low_level_arg=low_level_arg)

    @property
    def mapping_code(self):
        return f"{self.tracker_symbol.name}->getCounter()"

    @property
    def headers(self):
        return {'"lbm/inplace_streaming/TimestepTracker.h"'}


def generate_alternating_lbm_sweep(generation_context, class_name, collision_rule,
                                   lbm_config, lbm_optimisation=None,
                                   namespace='lbm', field_swaps=(), varying_parameters=(),
                                   inner_outer_split=False, ghost_layers_to_include=0,
                                   target=Target.CPU, data_type=None,
                                   cpu_openmp=None, cpu_vectorize_info=None, max_threads=None,
                                   **kernel_parameters):
    """Generates an Alternating lattice Boltzmann sweep class. This is in particular meant for
    in-place streaming patterns, but can of course also be used with two-fields patterns (why make it
    simple if you can make it complicated?). From a collision rule, two kernel implementations are
    generated; one for even, and one for odd timesteps. At run time, the correct one is selected
    according to a time step counter. This counter can be managed by the `walberla::lbm::TimestepTracker`
    class.

    Args:
        generation_context: See documentation of `pystencils_walberla.generate_sweep`
        class_name: Name of the generated class
        collision_rule: LbCollisionRule as returned by `lbmpy.create_lb_collision_rule`.
        lbm_config: configuration of the LB method. See lbmpy.LBMConfig
        lbm_optimisation: configuration of the optimisations of the LB method. See lbmpy.LBMOptimisation
        namespace: see documentation of `generate_sweep`
        field_swaps: see documentation of `generate_sweep`
        varying_parameters: see documentation of `generate_sweep`
        inner_outer_split: see documentation of `generate_sweep`
        ghost_layers_to_include: see documentation of `generate_sweep`
        target: An pystencils Target to define cpu or gpu code generation. See pystencils.Target
        data_type: default datatype for the kernel creation. Default is double
        cpu_openmp: if loops should use openMP or not.
        cpu_vectorize_info: dictionary containing necessary information for the usage of a SIMD instruction set.
        max_threads: only relevant for GPU kernels. Will be argument of `__launch_bounds__`.
        kernel_parameters: other parameters passed to the creation of a pystencils.CreateKernelConfig
    """
    config = config_from_context(generation_context, target=target, data_type=data_type, cpu_openmp=cpu_openmp,
                                 cpu_vectorize_info=cpu_vectorize_info, **kernel_parameters)

    # Add the lbm collision rule to the config
    lbm_config = replace(lbm_config, collision_rule=collision_rule)
    even_lbm_config = replace(lbm_config, timestep=Timestep.EVEN)

    ast_even = create_lb_ast(lbm_config=even_lbm_config, lbm_optimisation=lbm_optimisation, config=config)
    ast_even.function_name = 'even'
    kernel_even = KernelCallNode(ast_even)

    if is_inplace(lbm_config.streaming_pattern):
        odd_lbm_config = replace(lbm_config, timestep=Timestep.ODD)

        ast_odd = create_lb_ast(lbm_config=odd_lbm_config, lbm_optimisation=lbm_optimisation, config=config)
        ast_odd.function_name = 'odd'
        kernel_odd = KernelCallNode(ast_odd)
    else:
        kernel_odd = kernel_even

    tree = EvenIntegerCondition('timestep', kernel_even, kernel_odd, np.uint8)
    interface_mappings = [TimestepTrackerMapping(tree.parameter_symbol)]

    vec_info = config.cpu_vectorize_info
    openmp = config.cpu_openmp

    generate_selective_sweep(generation_context, class_name, tree,
                             interface_mappings=interface_mappings,
                             target=target, namespace=namespace,
                             field_swaps=field_swaps, varying_parameters=varying_parameters,
                             inner_outer_split=inner_outer_split, ghost_layers_to_include=ghost_layers_to_include,
                             cpu_vectorize_info=vec_info, cpu_openmp=openmp, max_threads=max_threads)
