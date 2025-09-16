# This file is part of waLBerla. waLBerla is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# waLBerla is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

from dataclasses import dataclass

from pystencils import (
    Field,
    FieldType,
    Assignment,
    CreateKernelConfig,
    DynamicType,
    Target,
)
from pystencils.stencil import offset_to_direction_string
from pystencils.codegen.properties import FieldBasePtr
from lbmpy import Stencil, LBStencil

from pystencilssfg import SfgComposer
from pystencilssfg.composer.basic_composer import KernelsAdder
from pystencilssfg.composer.custom import CustomGenerator
from pystencilssfg.ir import SfgKernelHandle, SfgEmptyNode
from pystencilssfg.ir.call_tree import SfgCallTreeNode
from pystencilssfg.ir.postprocessing import PostProcessingContext, SfgDeferredNode
from pystencilssfg.lang import SfgKernelParamVar, AugExpr, strip_ptr_ref
from pystencilssfg.lang.cpp import std
from pystencilssfg.lang.gpu import CudaAPI, HipAPI, ProvidesGpuRuntimeAPI

from ..api import GpuFieldPtr, Direction, CellInterval, uint_t
from ..build_config import get_build_config


@dataclass
class PackingKernelsContext:
    sfg: SfgComposer
    kns: KernelsAdder
    cfg: CreateKernelConfig


@dataclass
class CaptureBufferPointer(SfgDeferredNode):
    sfg: SfgComposer
    buffer_name: str
    buffer_span: std.span

    def expand(self, ppc: PostProcessingContext) -> SfgCallTreeNode:
        for param in ppc.live_variables:
            if (
                isinstance(param, SfgKernelParamVar)
                and param.wrapped.fields[0].name == self.buffer_name
                and param.wrapped.get_properties(FieldBasePtr)
            ):
                return self.sfg.init(param)(
                    AugExpr.format("{}.data()", self.buffer_span)
                )

        return SfgEmptyNode()


class GpuPdfFieldPackInfo(CustomGenerator):
    """Pack Info for lattice Boltzmann Gpu PDF fields.

    Generate a ghost layer exchange pack info for communicating lattice Boltzmann populations
    streaming across a block boundary,
    for use with `gpu::GpuField` and `gpu::UniformGpuScheme`.

    For a given velocity set, this pack info will only communicate those populations f_i
    from block A to a neighbor block B which are being advected, by the streaming step,
    from a cell in A to an adjacent cell in B.

    .. note::
        For the time being, this pack info is restricted to the *pull* streaming pattern.

    Args:
        name: Name of the generated pack info class
        stencil: Velocity set of the lattice Boltzmann method
        field: Symbolic representation of the PDF field
    """

    def __init__(self, name: str, stencil: LBStencil, field: Field):
        if field.index_dimensions > 1:
            raise ValueError(
                "GpuFieldPackInfo currently does not support higher-order tensor fields"
            )

        if isinstance(field.dtype, DynamicType):
            raise ValueError(
                "Cannot generate GpuFieldPackInfo for a dynamically-typed field"
            )

        self._name = name
        self._stencil = stencil
        self._full_stencil = (
            LBStencil(Stencil.D3Q27)
            if self._stencil.D == 3
            else LBStencil(Stencil.D2Q9)
        )

        #   Map storing the set of communicated populations for each communication direction
        self._communication_sets: dict[tuple[int, int, int], list[int]] = dict()
        for comm_dir in self._full_stencil:
            if indices := self._get_streaming_indices(comm_dir):
                self._communication_sets[comm_dir] = indices

        self._field = field
        self._dtype = field.dtype
        self._src_field = Field.new_field_with_different_name(
            self._field, f"{self._field.name}_src"
        )
        self._dst_field = Field.new_field_with_different_name(
            self._field, f"{self._field.name}_dst"
        )

    def generate(self, sfg: SfgComposer) -> None:
        base_class = f"walberla::sweepgen::UniformGpuFieldPackInfoBase< {self._name } >"
        sfg.include(
            "walberla/sweepgen/communication/UniformGpuFieldPackInfoBase.hpp"
        )

        build_config = get_build_config(sfg)

        pkc = PackingKernelsContext(
            sfg,
            kns=sfg.kernel_namespace(f"{self._name}_kernels"),
            cfg=build_config.get_pystencils_config(),
        )

        GpuAPI: type[ProvidesGpuRuntimeAPI]
        match pkc.cfg.get_target():
            case Target.CUDA:
                GpuAPI = CudaAPI
            case Target.HIP:
                GpuAPI = HipAPI
            case other:
                raise ValueError(
                    f"Invalid target for generating GpuFieldPackInfo: {other}"
                )

        pack_kernels: dict[tuple[int, int, int], SfgKernelHandle] = dict()
        unpack_kernels: dict[tuple[int, int, int], SfgKernelHandle] = dict()
        local_copy_kernels: dict[tuple[int, int, int], SfgKernelHandle] = dict()

        comm_dirs = self._communication_sets.keys()

        for comm_dir in comm_dirs:
            if not all(c == 0 for c in comm_dir):
                pack_kernels[comm_dir] = self._do_pack(pkc, comm_dir)
                unpack_kernels[comm_dir] = self._do_unpack(pkc, comm_dir)
                local_copy_kernels[comm_dir] = self._do_local_copy(pkc, comm_dir)

        src_gpu_field = GpuFieldPtr.create(self._src_field)
        gpu_field_type = strip_ptr_ref(src_gpu_field.get_dtype())
        dst_gpu_field = GpuFieldPtr.create(self._dst_field)
        buffer_span = std.span(self._dtype).var("buffer")
        dir = Direction().var("dir")
        src_interval = CellInterval(const=True, ref=True).var("srcInterval")
        dst_interval = CellInterval(const=True, ref=True).var("dstInterval")

        stream = GpuAPI.stream_t().var("stream")

        sfg.klass(self._name, bases=[f"public {base_class}"])(
            sfg.public(
                f"using Base = {base_class};",
                "using Base::Base;",
                f"using Field_T = {gpu_field_type.c_string()};",
                sfg.method("doPack").params(
                    src_gpu_field, buffer_span, dir, src_interval, stream
                )(
                    sfg.map_field(
                        self._src_field, src_gpu_field.with_cell_interval(src_interval)
                    ),
                    CaptureBufferPointer(sfg, "buffer", buffer_span),
                    sfg.switch(dir)
                    .cases(
                        {
                            Direction.from_offset(comm_dir): sfg.gpu_invoke(
                                pack_kernels[comm_dir], stream=stream
                            )
                            for comm_dir in comm_dirs
                        }
                    )
                    .default("/* unreachable */"),
                ),
                sfg.method("doUnpack").params(
                    src_gpu_field, buffer_span, dir, dst_interval, stream
                )(
                    sfg.map_field(
                        self._dst_field, src_gpu_field.with_cell_interval(dst_interval)
                    ),
                    CaptureBufferPointer(sfg, "buffer", buffer_span),
                    sfg.switch(dir)
                    .cases(
                        {
                            Direction.from_offset(comm_dir): sfg.gpu_invoke(
                                unpack_kernels[comm_dir], stream=stream
                            )
                            for comm_dir in comm_dirs
                        }
                    )
                    .default("/* unreachable */"),
                ),
                sfg.method("doLocalCopy").params(
                    src_gpu_field,
                    src_interval,
                    dst_gpu_field,
                    dst_interval,
                    dir,
                    stream,
                )(
                    sfg.map_field(
                        self._src_field, src_gpu_field.with_cell_interval(src_interval)
                    ),
                    sfg.map_field(
                        self._dst_field, dst_gpu_field.with_cell_interval(dst_interval)
                    ),
                    sfg.switch(dir)
                    .cases(
                        {
                            f"walberla::stencil::Direction::{offset_to_direction_string(comm_dir)}": sfg.gpu_invoke(
                                local_copy_kernels[comm_dir], stream=stream
                            )
                            for comm_dir in comm_dirs
                        }
                    )
                    .default("/* unreachable */"),
                ),
                sfg.method("elementsPerCell")
                .inline()
                .const()
                .params(dir)
                .returns(uint_t)(
                    sfg.switch(dir, autobreak=False).cases({
                        Direction.from_offset(comm_dir): f"return {len(elems)};"
                        for comm_dir, elems in self._communication_sets.items()
                    }).default("return 0;")
                ),
            )
        )

    def _pack_accesses(self, comm_dir: tuple[int, int, int]) -> list[Field.Access]:
        return [self._src_field.center(i) for i in self._communication_sets[comm_dir]]

    def _unpack_accesses(self, comm_dir: tuple[int, int, int]) -> list[Field.Access]:
        return [self._dst_field.center(i) for i in self._communication_sets[comm_dir]]

    def _get_streaming_indices(self, comm_dir) -> list[int]:
        if all(d == 0 for d in comm_dir):
            return []
        else:
            from lbmpy.advanced_streaming.communication import _extend_dir

            directions = set(_extend_dir(comm_dir)) & set(self._stencil)
            indices = sorted(self._stencil.index(d) for d in directions)
            return indices

    def _do_pack(
        self, pkc: PackingKernelsContext, comm_dir: tuple[int, int, int]
    ) -> SfgKernelHandle:
        pack_accs = self._pack_accesses(comm_dir)
        buffer = self._buffer(len(pack_accs))
        asms = [Assignment(buffer(i), acc) for i, acc in enumerate(pack_accs)]
        dir_str = offset_to_direction_string(comm_dir)
        return pkc.kns.create(asms, f"pack{dir_str}", pkc.cfg)

    def _do_unpack(
        self, pkc: PackingKernelsContext, comm_dir: tuple[int, int, int]
    ) -> SfgKernelHandle:
        unpack_accs = self._unpack_accesses(comm_dir)
        buffer = self._buffer(len(unpack_accs))
        asms = [Assignment(acc, buffer(i)) for i, acc in enumerate(unpack_accs)]
        dir_str = offset_to_direction_string(comm_dir)
        return pkc.kns.create(asms, f"unpack{dir_str}", pkc.cfg)

    def _do_local_copy(
        self, pkc: PackingKernelsContext, comm_dir: tuple[int, int, int]
    ) -> SfgKernelHandle:
        pack_accs = self._pack_accesses(comm_dir)
        unpack_accs = self._unpack_accesses(comm_dir)

        asms = [Assignment(dst, src) for dst, src in zip(unpack_accs, pack_accs)]
        dir_str = offset_to_direction_string(comm_dir)
        return pkc.kns.create(asms, f"localCopy{dir_str}", pkc.cfg)

    def _buffer(self, num_elems: int):
        return Field.create_generic(
            "buffer",
            1,
            field_type=FieldType.BUFFER,
            dtype=self._field.dtype,
            index_shape=(num_elems,),
        )
