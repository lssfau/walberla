import numpy as np
from pystencils_walberla.codegen import generate_selective_sweep, get_vectorize_instruction_set
from pystencils_walberla.kernel_selection import (
    AbstractInterfaceArgumentMapping, AbstractConditionNode, KernelCallNode)
from pystencils import TypedSymbol
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
        self.high_level_args = (self.tracker_symbol,)
        self.low_level_arg = low_level_arg

    @property
    def mapping_code(self):
        return f"{self.tracker_symbol.name}->getCounter()"

    @property
    def headers(self):
        return {'"lbm/inplace_streaming/TimestepTracker.h"'}


def generate_alternating_lbm_sweep(generation_context, class_name, collision_rule, streaming_pattern,
                                   namespace='lbm', field_swaps=(), varying_parameters=(),
                                   inner_outer_split=False, ghost_layers_to_include=0, optimization=None,
                                   **create_ast_params):
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
        streaming_pattern: Name of the streaming pattern; see `lbmpy.advanced_streaming`
        namespace: see documentation of `generate_sweep`
        field_swaps: see documentation of `generate_sweep`
        varying_parameters: see documentation of `generate_sweep`
        inner_outer_split: see documentation of `generate_sweep`
        ghost_layers_to_include: see documentation of `generate_sweep`
        optimization: dictionary containing optimization parameters, c.f. `lbmpy.creationfunctions`
        create_ast_params: Further parameters passed to `create_lb_ast`
    """
    optimization = default_lbm_optimization_parameters(generation_context, optimization)

    target = optimization['target']
    if not generation_context.cuda and target == 'gpu':
        return

    ast_even = create_lb_ast(collision_rule=collision_rule, streaming_pattern=streaming_pattern,
                             timestep=Timestep.EVEN, optimization=optimization, **create_ast_params)
    ast_even.function_name = 'even'
    kernel_even = KernelCallNode(ast_even)

    if is_inplace(streaming_pattern):
        ast_odd = create_lb_ast(collision_rule=collision_rule, streaming_pattern=streaming_pattern,
                                timestep=Timestep.ODD, optimization=optimization, **create_ast_params)
        ast_odd.function_name = 'odd'
        kernel_odd = KernelCallNode(ast_odd)
    else:
        kernel_odd = kernel_even

    tree = EvenIntegerCondition('timestep', kernel_even, kernel_odd, np.uint8)
    interface_mappings = [TimestepTrackerMapping(tree.parameter_symbol)]

    vec_info = optimization['vectorization']
    openmp = optimization['openmp']

    generate_selective_sweep(generation_context, class_name, tree,
                             interface_mappings=interface_mappings,
                             target=target, namespace=namespace,
                             field_swaps=field_swaps, varying_parameters=varying_parameters,
                             inner_outer_split=inner_outer_split, ghost_layers_to_include=ghost_layers_to_include,
                             cpu_vectorize_info=vec_info, cpu_openmp=openmp)


# ---------------------------------- Internal --------------------------------------------------------------------------


def default_lbm_optimization_parameters(generation_context, params):
    if params is None:
        params = dict()

    params['target'] = params.get('target', 'cpu')
    params['double_precision'] = params.get('double_precision', generation_context.double_accuracy)
    params['openmp'] = params.get('cpu_openmp', generation_context.openmp)
    params['vectorization'] = params.get('vectorization', {})

    if isinstance(params['vectorization'], bool):
        do_vectorization = params['vectorization']
        params['vectorization'] = dict()
    else:
        do_vectorization = True

    vec = params['vectorization']
    if isinstance(vec, dict):
        default_vec_is = get_vectorize_instruction_set(generation_context) if do_vectorization else None

        vec['instruction_set'] = vec.get('instruction_set', default_vec_is)
        vec['assume_inner_stride_one'] = vec.get('assume_inner_stride_one', True)
        vec['assume_aligned'] = vec.get('assume_aligned', False)
        vec['nontemporal'] = vec.get('nontemporal', False)
    return params
