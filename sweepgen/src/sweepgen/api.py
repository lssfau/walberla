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

from typing import cast

from pystencils import Field, DynamicType, FieldType, Target
from pystencils.types import (
    UserTypeSpec,
    create_type,
    PsCustomType,
    PsPointerType,
    PsType,
    PsStructType,
)
from pystencilssfg.lang import (
    AugExpr,
    SupportsFieldExtraction,
    SupportsVectorExtraction,
    Ref,
    ExprLike,
    cpptype,
    CppClass,
)
from pystencilssfg.lang.types import CppTypeFactory, CppType
from pystencilssfg.lang.cpp import std


real_t = PsCustomType("walberla::real_t")
cell_idx_t = PsCustomType("walberla::cell_idx_t")
uint_t = PsCustomType("walberla::uint_t")


class _PlainCppClass(AugExpr):
    _type: CppTypeFactory

    def __init__(self, const: bool = False, ref: bool = False):
        dtype = self._type(const=const, ref=ref)
        super().__init__(dtype)


class Vector2(AugExpr, SupportsVectorExtraction):
    _template = cpptype("walberla::Vector2< {element_type} >", "core/math/Vector2.h")

    def __init__(
        self, element_type: UserTypeSpec, const: bool = False, ref: bool = False
    ):
        self._element_type = create_type(element_type)
        dtype = self._template(element_type=element_type, const=const, ref=ref)
        super().__init__(dtype)

    def _extract_component(self, coordinate: int) -> AugExpr:
        if coordinate > 1:
            raise ValueError(f"Cannot extract component {coordinate} from Vector2")

        return AugExpr(self._element_type).bind("{}[{}]", self, coordinate)


class Vector3(AugExpr, SupportsVectorExtraction):
    _template = cpptype("walberla::Vector3< {element_type} >", "core/math/Vector3.h")

    def __init__(
        self, element_type: UserTypeSpec, const: bool = False, ref: bool = False
    ):
        self._element_type = create_type(element_type)
        dtype = self._template(element_type=element_type, const=const, ref=ref)
        super().__init__(dtype)

    def _extract_component(self, coordinate: int) -> AugExpr:
        if coordinate > 2:
            raise ValueError(f"Cannot extract component {coordinate} from Vector3")

        return AugExpr(self._element_type).bind("{}[{}]", self, coordinate)

    def __getitem__(self, idx: int | ExprLike):
        return AugExpr(self._element_type).bind("{}[{}]", self, idx)


class AABB(_PlainCppClass):
    _type = cpptype("walberla::AABB", "core/math/AABB.h")

    def min(self) -> Vector3:
        return Vector3(real_t).bind("{}.min()", self)

    def max(self) -> Vector3:
        return Vector3(real_t).bind("{}.max()", self)


class Cell(_PlainCppClass):
    _type = cpptype("walberla::Cell", "core/cell/Cell.h")

    def x(self) -> AugExpr:
        return AugExpr.format("{}.x()", self)

    def y(self) -> AugExpr:
        return AugExpr.format("{}.y()", self)

    def z(self) -> AugExpr:
        return AugExpr.format("{}.z()", self)

    def __getitem__(self, coord: int | ExprLike) -> AugExpr:
        return AugExpr.format("{}[{}]", self, coord)


class CellInterval(_PlainCppClass):
    _type = cpptype("walberla::CellInterval", "core/cell/CellInterval.h")

    def min(self) -> Cell:
        return Cell(ref=True).bind("{}.min()", self)

    def max(self) -> Cell:
        return Cell(ref=True).bind("{}.max()", self)


class Direction(_PlainCppClass):
    _type = cpptype("walberla::stencil::Direction", "stencil/Directions.h")

    @staticmethod
    def from_offset(offset: tuple[int, int, int]) -> str:
        from pystencils.stencil import offset_to_direction_string

        return f"walberla::stencil::Direction::{offset_to_direction_string(offset)}"


class BlockDataID(_PlainCppClass):
    _type = cpptype("walberla::BlockDataID", "domain_decomposition/BlockDataID.h")


class IBlock(_PlainCppClass):
    _type = cpptype("walberla::IBlock", "domain_decomposition/IBlock.h")

    def getData(self, dtype: str | PsType, id: BlockDataID) -> AugExpr:
        return AugExpr.format("{}.template getData< {} >({})", self, dtype, id)

    def getAABB(self) -> AABB:
        return AABB().bind("{}.getAABB()", self)


class IBlockPtr(_PlainCppClass):
    _type = cpptype("walberla::IBlock *", "domain_decomposition/IBlock.h")

    def getData(self, dtype: str | PsType, id: BlockDataID) -> AugExpr:
        return AugExpr.format("{}->template getData< {} >({})", self, dtype, id)

    def getAABB(self) -> AABB:
        return AABB().bind("{}->getAABB()", self)

    def deref(self) -> IBlock:
        return IBlock(ref=True).bind("*{}", self)


class SharedPtr(CppClass):
    template = cpptype("std::shared_ptr< {} >", "<memory>")


class StructuredBlockForest(AugExpr):
    typename = cpptype(
        "walberla::StructuredBlockForest", "blockforest/StructuredBlockForest.h"
    )

    def __init__(self, ref: bool = True, const: bool = False):
        dtype = self.typename(const=const, ref=ref)
        super().__init__(dtype)

    @staticmethod
    def shared_ptr(name: str = "blocks", const: bool = False):
        dtype = StructuredBlockForest.typename(const=const)
        return SharedPtr(dtype).var(name)

    @staticmethod
    def shared_ptr_ref(name: str = "blocks", const: bool = False):
        dtype = StructuredBlockForest.typename(const=const)
        return SharedPtr(dtype, const=True, ref=True).var(name)

    @staticmethod
    def ref(name: str = "blocks", const: bool = False):
        return StructuredBlockForest(ref=True).var(name)

    def getDomain(self) -> AABB:
        return AABB().bind("{}.getDomain()", self)

    def dx(self, level: int | AugExpr = 0) -> AugExpr:
        return AugExpr(real_t).bind("{}.dx({})", self, level)

    def dy(self, level: int | AugExpr = 0) -> AugExpr:
        return AugExpr(real_t).bind("{}.dy({})", self, level)

    def dz(self, level: int | AugExpr = 0) -> AugExpr:
        return AugExpr(real_t).bind("{}.dz({})", self, level)

    def cell_extents(self, coord: int, level: int | AugExpr = 0) -> AugExpr:
        dx_getter = [self.dx, self.dy, self.dz][coord]
        return dx_getter(level)

    def cellsPerBlock(self, coord: int, block: AugExpr) -> AugExpr:
        x = ["X", "Y", "Z"][coord]
        return AugExpr(uint_t).bind("{}.getNumberOf{}Cells({})", self, x, block)

    def getDomainCellBB(self, level: ExprLike) -> CellInterval:
        return CellInterval(ref=True).bind("{}.getDomainCellBB({})", self, level)

    def getBlockCellBB(self, block: IBlock) -> CellInterval:
        return CellInterval().bind("{}.getBlockCellBB({})", self, block)

    def getLevel(self, block: AugExpr) -> AugExpr:
        return AugExpr(uint_t).bind("{}.getLevel({})", self, block)


class GenericWalberlaField(AugExpr, SupportsFieldExtraction):
    """Common base class for GhostLayerField and GpuField defining their shared interface."""

    def __init__(
        self,
        element_type: PsType,
        field_type: PsType,
        spatial_rank: int,
        ptr: bool = False,
        ref: bool = False,
    ):
        if ptr and ref:
            raise ValueError("At most one of `ptr` and `ref` may be true.")

        self._element_type = element_type
        self._field_type = field_type
        self._spatial_rank = spatial_rank

        obj_type: PsType
        if ptr:
            obj_type = PsPointerType(self._field_type)
        elif ref:
            obj_type = Ref(self._field_type)
        else:
            obj_type = self._field_type

        super().__init__(obj_type)

        self._extraction = GhostLayerFieldExtraction(self, None)

    @property
    def _a(self) -> AugExpr:
        """Member access"""
        if isinstance(self.dtype, PsPointerType):
            return AugExpr.format("{}->", self)
        else:
            return AugExpr.format("{}.", self)

    @property
    def field_type(self) -> PsType:
        return self._field_type

    def _extract_ptr(self) -> AugExpr:
        return self._extraction._extract_ptr()

    def _extract_size(self, coordinate: int) -> AugExpr | None:
        return self._extraction._extract_size(coordinate)

    def _extract_stride(self, coordinate: int) -> AugExpr | None:
        return self._extraction._extract_stride(coordinate)

    def with_cell_interval(self, ci: CellInterval | None) -> GhostLayerFieldExtraction:
        return GhostLayerFieldExtraction(self, ci)

    def cloneUninitialized(self) -> AugExpr:
        return AugExpr.format("{}cloneUninitialized()", self._a)

    def swapDataPointers(self, other: AugExpr) -> AugExpr:
        return AugExpr.format("{}swapDataPointers({});", self._a, other)


class GhostLayerFieldPtr(GenericWalberlaField):
    _template = cpptype(
        "walberla::field::GhostLayerField< {element_type}, {fsize} >",
        "field/GhostLayerField.h",
    )

    @staticmethod
    def create(field: Field, const: bool = False):
        if field.index_dimensions > 1:
            raise ValueError(
                "Cannot map fields with more than one index dimension to field::GhostLayerField."
            )

        element_type = field.dtype

        if isinstance(element_type, DynamicType):
            raise ValueError(
                f"Cannot create GhostLayerField from dynamically typed field {field}"
            )

        fsize = field.index_shape[0] if field.index_shape else 1

        return GhostLayerFieldPtr(
            element_type, field.spatial_dimensions, fsize, const=const
        ).var(field.name)

    def __init__(
        self,
        element_type: UserTypeSpec,
        spatial_rank: int,
        fsize: int,
        const: bool = False,
    ):
        element_type = create_type(element_type)
        field_type = self._template(element_type=element_type, fsize=fsize, const=const)

        super().__init__(element_type, field_type, spatial_rank, ptr=True)


class GpuFieldPtr(GenericWalberlaField):
    _template = cpptype(
        "walberla::gpu::GPUField< {element_type} >",
        "gpu/GPUField.h",
    )

    @staticmethod
    def create(field: Field, const: bool = False):
        if field.index_dimensions > 1:
            raise ValueError(
                "Cannot map fields with more than one index dimension to gpu::GpuField."
            )

        element_type = field.dtype

        if isinstance(element_type, DynamicType):
            raise ValueError(
                f"Cannot create GpuField from dynamically typed field {field}"
            )

        fsize = field.index_shape[0] if field.index_shape else 1

        return GpuFieldPtr(
            element_type, field.spatial_dimensions, fsize, const=const
        ).var(field.name)

    def __init__(
        self,
        element_type: UserTypeSpec,
        spatial_rank: int,
        fsize: int,
        const: bool = False,
    ):
        element_type = create_type(element_type)
        field_type = self._template(element_type=element_type, const=const)

        super().__init__(element_type, field_type, spatial_rank, ptr=True)


class GhostLayerFieldExtraction(SupportsFieldExtraction):
    def __init__(
        self,
        field_ptr: GenericWalberlaField,
        cell_interval: AugExpr | None = None,
    ) -> None:
        self._field_ptr = field_ptr
        self._ci = cell_interval

        self._size_calls = ("xSize()", "ySize()", "zSize()")[
            : field_ptr._spatial_rank
        ] + ("fSize()",)
        self._stride_calls = ("xStride()", "yStride()", "zStride()")[
            : field_ptr._spatial_rank
        ] + ("fStride()",)

    def _extract_ptr(self) -> AugExpr:
        data_at: AugExpr | str
        if self._ci is not None:
            ci = self._ci
            data_at = AugExpr.format("{ci}.xMin(), {ci}.yMin(), {ci}.zMin(), 0", ci=ci)
        else:
            data_at = "0, 0, 0, 0"

        return AugExpr.format("{}->dataAt({})", self._field_ptr, data_at)

    def _extract_size(self, coordinate: int) -> AugExpr | None:
        if coordinate > self._field_ptr._spatial_rank:
            d = self._field_ptr._spatial_rank
            raise IndexError(
                f"Invalid coordinate for {d}-dimensional waLBerla field: {coordinate}"
            )

        size_call = self._size_calls[coordinate]

        if self._ci is not None and coordinate < 3:
            return AugExpr.format("{}.{}", self._ci, size_call)
        else:
            return AugExpr.format("{}->{}", self._field_ptr, size_call)

    def _extract_stride(self, coordinate: int) -> AugExpr | None:
        if coordinate > self._field_ptr._spatial_rank:
            d = self._field_ptr._spatial_rank
            raise IndexError(
                f"Invalid coordinate for {d}-dimensional waLBerla field: {coordinate}"
            )

        stride_call = self._stride_calls[coordinate]
        return AugExpr.format("{}->{}", self._field_ptr, stride_call)


def glfield(field: Field, ci: str | AugExpr | None = None):
    field_ptr = GhostLayerFieldPtr.create(field)

    if ci is not None and not isinstance(ci, AugExpr):
        ci = AugExpr(PsCustomType("walberla::CellInterval", const=True)).var(ci)

    assert not isinstance(ci, str)

    return GhostLayerFieldExtraction(field_ptr, ci)


class MemTags:
    _header = "walberla/experimental/Memory.hpp"

    host = cpptype("walberla::experimental::memory::memtag::host", _header)()
    unified = cpptype("walberla::experimental::memory::memtag::unified", _header)()


class SparseIndexList(AugExpr):
    _template = cpptype(
        "walberla::experimental::sweep::SparseIndexList< {IndexStruct}, {memtag_t} >",
        "walberla/experimental/sweep/SparseIndexList.hpp",
    )

    def __init__(
        self,
        idx_struct: PsStructType,
        memtag_t: PsType = MemTags.host,
        const: bool = False,
        ref: bool = False,
    ):
        self._idx_struct = idx_struct
        dtype = self._template(
            IndexStruct=idx_struct, memtag_t=memtag_t, const=const, ref=ref
        )
        super().__init__(dtype)

    def view_type(self) -> std.span:
        return std.span(self._idx_struct)

    def view(self, block: IBlock) -> std.span:
        return std.span(self._idx_struct).bind("{}.view({})", self, block)

    def bufferId(self) -> BlockDataID:
        return BlockDataID().bind("{}.bufferId()", self)

    @staticmethod
    def from_field(
        field: Field,
        target: Target | None = None,
        const: bool = False,
        ref: bool = False,
    ):
        if (
            not isinstance(field.dtype, PsStructType)
            or field.field_type != FieldType.INDEXED
        ):
            raise ValueError(
                "Can only create instances of SparseIndexList from index fields holding structures."
            )

        if target is not None and target.is_gpu():
            memtag = MemTags.unified
        else:
            memtag = MemTags.host

        return SparseIndexList(field.dtype, memtag_t=memtag, const=const, ref=ref).var(
            field.name
        )


class IndexListBufferPtr(AugExpr, SupportsFieldExtraction):
    _template = cpptype(
        "walberla::experimental::sweep::internal::IndexListBuffer< {IndexStruct} >",
        "walberla/experimental/sweep/SparseIndexList.hpp",
    )

    def __init__(self, idx_struct: PsStructType):
        self._idx_struct = idx_struct
        dtype = self._template(IndexStruct=idx_struct)
        self._base_type = dtype
        super().__init__(PsPointerType(dtype))

    @property
    def field_type(self) -> CppType:
        return cast(CppType, self._base_type)

    def pointerCpu(self) -> AugExpr:
        return AugExpr().format("{}->pointerCpu()", self)

    def vector(self) -> AugExpr:
        return AugExpr().format("{}->vector()", self)

    def _extract_ptr(self) -> AugExpr:
        return self.pointerCpu()

    def _extract_size(self, coordinate: int) -> AugExpr | None:
        if coordinate == 0:
            return AugExpr.format("{}.size()", self.vector())
        elif coordinate == 1:
            return AugExpr.format("1")
        else:
            raise ValueError()

    def _extract_stride(self, coordinate: int) -> AugExpr | None:
        if coordinate <= 1:
            return AugExpr.format("1")
        else:
            raise ValueError()


CellIdx = PsStructType(
    (
        ("x", create_type("int64_t")),
        ("y", create_type("int64_t")),
        ("z", create_type("int64_t")),
    ),
    "walberla::experimental::sweep::CellIdx",
)
