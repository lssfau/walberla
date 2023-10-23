from dataclasses import replace
from itertools import product

import numpy as np
import sympy as sp

from jinja2 import Environment, PackageLoader, StrictUndefined

from pystencils import Assignment, CreateKernelConfig, create_kernel, Field, FieldType, fields, Target
from pystencils.stencil import offset_to_direction_string
from pystencils.typing import TypedSymbol
from pystencils.stencil import inverse_direction
from pystencils.bit_masks import flag_cond

from lbmpy.advanced_streaming import get_accessor, is_inplace, get_timesteps, Timestep
from lbmpy.advanced_streaming.communication import _extend_dir
from lbmpy.enums import Stencil
from lbmpy.stencils import LBStencil

from pystencils_walberla.cmake_integration import CodeGenerationContext
from pystencils_walberla.kernel_selection import KernelFamily, KernelCallNode, SwitchNode
from pystencils_walberla.jinja_filters import add_pystencils_filters_to_jinja_env
from pystencils_walberla.utility import config_from_context

from lbmpy_walberla.alternating_sweeps import EvenIntegerCondition
from lbmpy_walberla.utility import timestep_suffix


def generate_packing_kernels(generation_context: CodeGenerationContext, class_name: str,
                             stencil: LBStencil, streaming_pattern: str = 'pull',
                             namespace='lbm', nonuniform: bool = False,
                             target: Target = Target.CPU, data_type=None, cpu_openmp: bool = False,
                             **create_kernel_params):

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

    jinja_context = {
        'class_name': class_name,
        'namespace': namespace,
        'nonuniform': nonuniform,
        'target': target.name.lower(),
        'dtype': "float" if is_float else "double",
        'is_gpu': target == Target.GPU,
        'kernels': kernels,
        'inplace': is_inplace(streaming_pattern),
        'direction_sizes': cg.get_direction_sizes(),
        'stencil_size': stencil.Q,
        'dimension': stencil.D,
        'src_field': cg.src_field,
        'dst_field': cg.dst_field
    }

    if nonuniform:
        jinja_context['mask_field'] = cg.mask_field

    template_name = "NonuniformPackingKernels" if nonuniform else "PackingKernels"

    env = Environment(loader=PackageLoader('lbmpy_walberla'), undefined=StrictUndefined)
    add_pystencils_filters_to_jinja_env(env)
    header = env.get_template(f"{template_name}.tmpl.h").render(**jinja_context)
    source = env.get_template(f"{template_name}.tmpl.cpp").render(**jinja_context)

    source_extension = "cu" if target == Target.GPU and generation_context.cuda else "cpp"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)


#   ------------------------------ INTERNAL ----------------------------------------------------------------------------

class PackingKernelsCodegen:

    def __init__(self, stencil, streaming_pattern, class_name, config: CreateKernelConfig,
                 src_field=None, dst_field=None):
        self.stencil = stencil
        self.dim = stencil.D
        self.values_per_cell = stencil.Q
        self.full_stencil = LBStencil(Stencil.D3Q27) if self.dim == 3 else LBStencil(Stencil.D2Q9)
        self.streaming_pattern = streaming_pattern
        self.inplace = is_inplace(streaming_pattern)
        self.class_name = class_name
        self.config = config
        self.data_type = config.data_type['pdfs'].numpy_dtype

        self.src_field = src_field if src_field else fields(f'pdfs_src({stencil.Q}) :{self.data_type}[{stencil.D}D]')
        self.dst_field = dst_field if dst_field else fields(f'pdfs_dst({stencil.Q}) :{self.data_type}[{stencil.D}D]')

        self.accessors = [get_accessor(streaming_pattern, t) for t in get_timesteps(streaming_pattern)]
        self.mask_field = fields(f'mask : uint32 [{self.dim}D]', layout=src_field.layout)

    def create_uniform_kernel_families(self, kernels_dict=None):
        kernels = dict() if kernels_dict is None else kernels_dict

        kernels['packAll'] = self.get_pack_all_kernel_family()
        kernels['unpackAll'] = self.get_unpack_all_kernel_family()
        kernels['localCopyAll'] = self.get_local_copy_all_kernel_family()

        kernels['packDirection'] = self.get_pack_direction_kernel_family()
        kernels['unpackDirection'] = self.get_unpack_direction_kernel_family()
        kernels['localCopyDirection'] = self.get_local_copy_direction_kernel_family()
        return kernels

    def create_nonuniform_kernel_families(self, kernels_dict=None):
        kernels = dict() if kernels_dict is None else kernels_dict
        kernels['unpackRedistribute'] = self.get_unpack_redistribute_kernel_family()
        kernels['packPartialCoalescence'] = self.get_pack_partial_coalescence_kernel_family()
        kernels['zeroCoalescenceRegion'] = self.get_zero_coalescence_region_kernel_family()
        kernels['unpackCoalescence'] = self.get_unpack_coalescence_kernel_family()

        return kernels

    # --------------------------- Pack / Unpack / LocalCopy All --------------------------------------------------------

    def get_pack_all_ast(self, timestep):
        config = replace(self.config, ghost_layers=0)

        buffer = self._buffer(self.values_per_cell)
        src, _ = self._stream_out_accs(timestep)
        assignments = [Assignment(buffer(i), src[i]) for i in range(self.values_per_cell)]
        ast = create_kernel(assignments, config=config)
        ast.function_name = 'pack_ALL' + timestep_suffix(timestep)
        return ast

    def get_pack_all_kernel_family(self):
        if not self.inplace:
            tree = KernelCallNode(self.get_pack_all_ast(Timestep.BOTH))
        else:
            even_call = KernelCallNode(self.get_pack_all_ast(Timestep.EVEN))
            odd_call = KernelCallNode(self.get_pack_all_ast(Timestep.ODD))
            tree = EvenIntegerCondition('timestep', even_call, odd_call, parameter_dtype=np.uint8)
        return KernelFamily(tree, self.class_name)

    def get_unpack_all_ast(self, timestep):
        config = replace(self.config, ghost_layers=0)

        buffer = self._buffer(self.values_per_cell)
        _, dst = self._stream_out_accs(timestep)
        assignments = [Assignment(dst[i], buffer(i)) for i in range(self.values_per_cell)]
        ast = create_kernel(assignments, config=config)
        ast.function_name = 'unpack_ALL' + timestep_suffix(timestep)
        return ast

    def get_unpack_all_kernel_family(self):
        if not self.inplace:
            tree = KernelCallNode(self.get_unpack_all_ast(Timestep.BOTH))
        else:
            even_call = KernelCallNode(self.get_unpack_all_ast(Timestep.EVEN))
            odd_call = KernelCallNode(self.get_unpack_all_ast(Timestep.ODD))
            tree = EvenIntegerCondition('timestep', even_call, odd_call, parameter_dtype=np.uint8)
        return KernelFamily(tree, self.class_name)

    def get_local_copy_all_ast(self, timestep):
        config = replace(self.config, ghost_layers=0)

        src, dst = self._stream_out_accs(timestep)
        assignments = [Assignment(dst[i], src[i]) for i in range(self.values_per_cell)]
        ast = create_kernel(assignments, config=config)
        ast.function_name = 'localCopy_ALL' + timestep_suffix(timestep)
        return ast

    def get_local_copy_all_kernel_family(self):
        if not self.inplace:
            tree = KernelCallNode(self.get_local_copy_all_ast(Timestep.BOTH))
        else:
            even_call = KernelCallNode(self.get_local_copy_all_ast(Timestep.EVEN))
            odd_call = KernelCallNode(self.get_local_copy_all_ast(Timestep.ODD))
            tree = EvenIntegerCondition('timestep', even_call, odd_call, parameter_dtype=np.uint8)
        return KernelFamily(tree, self.class_name)

    # --------------------------- Pack / Unpack / LocalCopy Direction --------------------------------------------------

    def get_pack_direction_ast(self, comm_dir, timestep):
        config = replace(self.config, ghost_layers=0)

        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(comm_dir)
        buffer = self._buffer(len(streaming_dirs))
        src, _ = self._stream_out_accs(timestep)
        assignments = []
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        for i, d in enumerate(dir_indices):
            assignments.append(Assignment(buffer(i), src[d]))
        ast = create_kernel(assignments, config=config)
        ast.function_name = f'pack_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_pack_direction_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_pack_direction_ast)

    def get_unpack_direction_ast(self, comm_dir, timestep):
        config = replace(self.config, ghost_layers=0)

        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(inverse_direction(comm_dir))
        buffer = self._buffer(len(streaming_dirs))
        _, dst = self._stream_out_accs(timestep)
        assignments = []
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        for i, d in enumerate(dir_indices):
            assignments.append(Assignment(dst[d], buffer(i)))
        ast = create_kernel(assignments, config=config)
        ast.function_name = f'unpack_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_unpack_direction_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_unpack_direction_ast)

    def get_local_copy_direction_ast(self, comm_dir, timestep):
        config = replace(self.config, ghost_layers=0)

        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(comm_dir)
        src, dst = self._stream_out_accs(timestep)
        assignments = []
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        for direction in dir_indices:
            assignments.append(Assignment(dst[direction], src[direction]))
        ast = create_kernel(assignments, config=config)
        ast.function_name = f'localCopy_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_local_copy_direction_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_local_copy_direction_ast)

    # --------------------------- Pack / Unpack / LocalCopy Coarse to Fine ---------------------------------------------

    def get_unpack_redistribute_ast(self, comm_dir, timestep):
        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(inverse_direction(comm_dir))
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        buffer = self._buffer(self.values_per_cell)
        _, dst = self._stream_out_accs(timestep)
        orthos = self.orthogonal_principals(comm_dir)
        sub_dirs = self.contained_principals(comm_dir)
        orthogonal_combinations = self.linear_combinations(orthos)
        subdir_combinations = self.linear_combinations_nozero(sub_dirs)
        second_gl_dirs = [o + s for o, s in product(orthogonal_combinations, subdir_combinations)]
        negative_dir_correction = np.array([(1 if d == -1 else 0) for d in comm_dir])
        assignments = []
        for offset in orthogonal_combinations:
            o = offset + negative_dir_correction
            for d in range(self.values_per_cell):
                field_acc = dst[d].get_shifted(*o)
                assignments.append(Assignment(field_acc, buffer(d)))

        for offset in second_gl_dirs:
            o = offset + negative_dir_correction
            for d in dir_indices:
                field_acc = dst[d].get_shifted(*o)
                assignments.append(Assignment(field_acc, buffer(d)))

        function_name = f'unpackRedistribute_{dir_string}' + timestep_suffix(timestep)
        iteration_slice = tuple(slice(None, None, 2) for _ in range(self.dim))
        config = CreateKernelConfig(function_name=function_name, iteration_slice=iteration_slice,
                                    data_type=self.data_type, ghost_layers=0, allow_double_writes=True,
                                    cpu_openmp=self.config.cpu_openmp, target=self.config.target)

        return create_kernel(assignments, config=config)

    def get_unpack_redistribute_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_unpack_redistribute_ast)

    def get_local_copy_redistribute_ast(self, comm_dir, timestep):
        #   TODO
        raise NotImplementedError()

    def get_local_copy_redistribute_kernel_family(self):
        #   TODO
        raise NotImplementedError()

    # --------------------------- Pack / Unpack / LocalCopy Fine to Coarse ---------------------------------------------

    def get_pack_partial_coalescence_ast(self, comm_dir, timestep):
        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(comm_dir)
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        buffer = self._buffer(self.values_per_cell)
        src, _ = self._stream_in_accs(timestep.next())
        mask = self.mask_field

        offsets = list(product(*((0, 1) for _ in comm_dir)))
        assignments = []
        for i, d in enumerate(dir_indices):
            acc = 0
            for o in offsets:
                acc += flag_cond(d, mask[o], src[d].get_shifted(*o))
            assignments.append(Assignment(buffer(i), acc))

        iteration_slice = tuple(slice(None, None, 2) for _ in range(self.dim))
        config = replace(self.config, iteration_slice=iteration_slice, ghost_layers=0)

        ast = create_kernel(assignments, config=config)
        ast.function_name = f'packPartialCoalescence_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_pack_partial_coalescence_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_pack_partial_coalescence_ast)

    def get_unpack_coalescence_ast(self, comm_dir, timestep):
        config = replace(self.config, ghost_layers=0)

        assert not all(d == 0 for d in comm_dir)
        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(inverse_direction(comm_dir))
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        buffer = self._buffer(self.values_per_cell)
        _, dst = self._stream_in_accs(timestep.next())

        coalescence_factor = sp.Rational(1, 2 ** self.dim)

        assignments = []
        for i, d in enumerate(dir_indices):
            assignments.append(Assignment(dst[d], dst[d] + coalescence_factor * buffer(i)))

        ast = create_kernel(assignments, config=config)
        ast.function_name = f'unpackCoalescence_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_unpack_coalescence_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_unpack_coalescence_ast)

    def get_zero_coalescence_region_ast(self, comm_dir, timestep):
        config = replace(self.config, ghost_layers=0)

        dir_string = offset_to_direction_string(comm_dir)
        streaming_dirs = self.get_streaming_dirs(inverse_direction(comm_dir))
        dir_indices = sorted(self.stencil.index(d) for d in streaming_dirs)
        if len(dir_indices) == 0:
            return None
        _, dst = self._stream_in_accs(timestep.next())

        assignments = []
        for i, d in enumerate(dir_indices):
            assignments.append(Assignment(dst[d], 0.0))

        ast = create_kernel(assignments, config=config)
        ast.function_name = f'zeroCoalescenceRegion_{dir_string}' + timestep_suffix(timestep)
        return ast

    def get_zero_coalescence_region_kernel_family(self):
        return self._construct_directionwise_kernel_family(self.get_zero_coalescence_region_ast)

    #   TODO
    def get_local_copy_partial_coalescence_ast(self, comm_dir, timestep):
        raise NotImplementedError()

    def get_local_copy_partial_coalescence_kernel_family(self):
        raise NotImplementedError()

    # ------------------------------------------ Utility ---------------------------------------------------------------

    def get_streaming_dirs(self, comm_dir):
        if all(d == 0 for d in comm_dir):
            return set()
        else:
            return set(_extend_dir(comm_dir)) & set(self.stencil)

    def get_direction_sizes(self):
        return [len(self.get_streaming_dirs(d)) for d in self.full_stencil]

    def principal(self, i):
        e_i = np.zeros(self.dim, dtype=int)
        e_i[i] = 1
        return e_i

    def principals(self):
        """Returns the principal directions for the given dimension"""
        return tuple(self.principal(i) for i in range(self.dim))

    def orthogonal_principals(self, comm_dir):
        """Returns the positive principal directions orthogonal to the comm_dir"""
        return tuple(p for i, p in enumerate(self.principals()) if comm_dir[i] == 0)

    def contained_principals(self, comm_dir):
        """Returns the (positive or negative) principal directions contained in comm_dir"""
        vecs = []
        for i, d in enumerate(comm_dir):
            if d != 0:
                vecs.append(d * self.principal(i))
        return vecs

    def linear_combinations(self, vectors):
        if not vectors:
            return [np.zeros(self.dim, dtype=int)]
        else:
            rest = self.linear_combinations(vectors[1:])
            return rest + [vectors[0] + r for r in rest]

    def linear_combinations_nozero(self, vectors):
        if len(vectors) == 1:
            return [vectors[0]]
        else:
            rest = self.linear_combinations_nozero(vectors[1:])
            return rest + [vectors[0]] + [vectors[0] + r for r in rest]

    # --------------------------- Private Members ----------------------------------------------------------------------

    def _construct_directionwise_kernel_family(self, create_ast_callback):
        subtrees = []
        direction_symbol = TypedSymbol('dir', dtype='stencil::Direction')
        for t in get_timesteps(self.streaming_pattern):
            cases_dict = dict()
            for comm_dir in self.full_stencil:
                if all(d == 0 for d in comm_dir):
                    continue
                dir_string = offset_to_direction_string(comm_dir)
                ast = create_ast_callback(comm_dir, t)
                if ast is None:
                    continue
                kernel_call = KernelCallNode(ast)
                cases_dict[f"stencil::{dir_string}"] = kernel_call
            subtrees.append(SwitchNode(direction_symbol, cases_dict))

        if not self.inplace:
            tree = subtrees[0]
        else:
            tree = EvenIntegerCondition('timestep', subtrees[Timestep.EVEN.idx], subtrees[Timestep.ODD.idx],
                                        parameter_dtype=np.uint8)
        return KernelFamily(tree, self.class_name)

    def _stream_out_accs(self, timestep):
        accessor = self.accessors[timestep.idx]
        src_stream_out_accs = accessor.write(self.src_field, self.stencil)
        dst_stream_out_accs = accessor.write(self.dst_field, self.stencil)
        return src_stream_out_accs, dst_stream_out_accs

    def _stream_in_accs(self, timestep):
        accessor = self.accessors[timestep.idx]
        src_stream_in_accs = accessor.read(self.src_field, self.stencil)
        dst_stream_in_accs = accessor.read(self.dst_field, self.stencil)
        return src_stream_in_accs, dst_stream_in_accs

    def _buffer(self, size):
        return Field.create_generic('buffer', spatial_dimensions=1, field_type=FieldType.BUFFER,
                                    dtype=self.data_type,
                                    index_shape=(size,))
