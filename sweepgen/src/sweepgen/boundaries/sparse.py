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

from abc import abstractmethod

from pystencils import (
    Field,
    FieldType,
    TypedSymbol,
    Assignment,
    AssignmentCollection,
    create_type,
)
from pystencils.types import PsStructType, PsCustomType

from lbmpy.methods import AbstractLbMethod
from lbmpy.boundaries.boundaryconditions import LbBoundary

from pystencilssfg import SfgComposer
from pystencilssfg.lang import HeaderFile, SfgVar
from pystencilssfg.composer.custom import CustomGenerator

from ..sweep import Sweep
from ..api import SparseIndexList, MemTags
from ..build_config import get_build_config
from ..core.properties import PropertiesContainerBuilder

BoundaryIndexType = create_type("int32")


class SparseBoundaryDefinition:

    def __init__(self, name: str) -> None:
        self._name = name

    @property
    def name(self) -> str:
        return self._name

    @abstractmethod
    def index_struct_type(self) -> PsStructType: ...

    @abstractmethod
    def factory_crtp(self) -> str: ...

    @abstractmethod
    def includes(self) -> set[HeaderFile]: ...

    @abstractmethod
    def __call__(
        self, f_out, f_in, dir_symbol, inv_dir, lb_method, index_field, force_vector
    ) -> Assignment | list[Assignment] | AssignmentCollection: ...


class SparseBoundary(CustomGenerator):
    """Generate a sparse boundary handling sweep from a given boundary definition.

    Generate a sparse boundary handling sweep that applies a given boundary condition
    on cells or lattice links listed on a given index vector.
    This generator supports boundary conditions provided by ``lbmpy.boundaries``
    as well as SweepGen-specific sparse boundaries deriving from `SparseBoundaryDefinition`.

    Sparse boundary sweeps rely on an index list to iterate boundary links efficiently.
    This index list is created by a factory class.
    The factory comprises two parts:
    A CRTP base class located in a header file of the SweepGen runtime library, which implements the actual
    index list construction logic and incorporates implementation details via templates;
    and a minimal, boundary-specific subclass that is generated alongside the boundary sweep.

    .. admonition:: Example Apps
        :class: seealso

        - :walberla-example-app:`ParallelPlates`
        - :walberla-example-app:`FlowAroundSphere`
    """

    def __init__(
        self,
        bc: LbBoundary | SparseBoundaryDefinition,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
    ):
        if isinstance(bc, LbBoundary):
            bc = LbmpyBoundaryWrapper(bc)

        self._bc: SparseBoundaryDefinition = bc
        self._name = self._bc.name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._stencil = lb_method.stencil

    def generate(self, sfg: SfgComposer) -> None:
        bc = self._bc

        for header in bc.includes():
            sfg.include(header)

        #   Get assignments for bc
        bc_asm = self._get_kernel_assignments()

        #   Build generator config
        bc_cfg = get_build_config(sfg).get_pystencils_config()
        bc_cfg.index_dtype = BoundaryIndexType
        index_field = self._index_field()
        bc_cfg.index_field = index_field

        if isinstance(bc, LbmpyBoundaryWrapper):
            bc.render_data_struct(sfg)
            bc.render_link_struct(sfg)

        #   Prepare sweep
        bc_sweep = Sweep(self._name, bc_asm, bc_cfg)

        #   Emit code
        sfg.generate(bc_sweep)

        #   Build factory
        factory_name = f"{self._name}Factory"
        factory_crtp_base = f"{bc.factory_crtp()}< {factory_name} >"
        memtag_t = MemTags.unified if bc_cfg.get_target().is_gpu() else MemTags.host

        index_vector_name = "indexVector"
        index_vector = SparseIndexList(
            bc.index_struct_type(), memtag_t=memtag_t, ref=True
        ).var(index_vector_name)

        #   Create a property cache for forwarding all sweep parameters
        #   except the index vector
        sweep_type = bc_sweep.generated_class()
        sweep_props_container = bc_sweep.properties_container

        factory_args_factory = PropertiesContainerBuilder()
        for param in sweep_props_container.ctor_params:
            if param.name != index_vector_name:
                factory_args_factory.add_property(param, getter=False, setter=False)

        with sfg.namespace("detail"):
            factory_args_struct = factory_args_factory.render_struct(
                sfg, f"{self._name}FactoryArgs"
            )
            factory_args_cache = factory_args_struct().var("factory_args_")
            factory_args_cache_ref = factory_args_struct().bind(
                f"this->{factory_args_cache}"
            )

        #   Factory constructor

        blockforest_ptr = PropertiesContainerBuilder.blockforest_shared_ptr
        args_without_blockforest = tuple(
            filter(lambda p: p != blockforest_ptr, factory_args_struct.ctor_params)
        )
        factory_ctor = (
            sfg.constructor(blockforest_ptr, *args_without_blockforest)
            .init("Base")(blockforest_ptr)
            .init(factory_args_cache)(*factory_args_cache.ctor_params)
        )

        sweep_ctor_args = {
            index_field.name: index_vector,
        } | {
            p.name: factory_args_cache_ref.get(p.var)
            for p in factory_args_struct.properties.values()
        }

        stencil_name = self._lb_method.stencil.name
        sfg.include(f"stencil/{stencil_name}.h")

        sfg.klass(factory_name, bases=[f"public {factory_crtp_base}"])(
            sfg.public(
                f"using Base = {factory_crtp_base};",
                f"friend class {factory_crtp_base};",
                f"using Stencil = walberla::stencil::{stencil_name};",
                f"using Sweep = {sweep_type.get_dtype().c_string()};",
                f"using memtag_t = {memtag_t.c_string()};",
                f"using IdxStruct = {bc.index_struct_type().name};",
                (
                    f"using DataStruct = {bc.data_struct_name()};"
                    if isinstance(bc, LbmpyBoundaryWrapper) and bc.generate_add_data()
                    else "using DataStruct =  void;"
                ),
                factory_ctor,
            ),
            sfg.private(
                factory_args_cache,
                sfg.method(
                    "irregularFromIndexVector",
                )
                .returns(sweep_type.get_dtype())
                .inline()(sfg.expr("return {};", sweep_type.ctor(**sweep_ctor_args))),
            ),
        )

    def _index_field(self):
        return Field(
            "indexVector",
            FieldType.INDEXED,
            self._bc.index_struct_type(),
            (0,),
            (TypedSymbol("indexVectorLength", BoundaryIndexType), 1),
            (1, 1),
        )

    def _get_kernel_assignments(self) -> list[Assignment]:
        index_field = self._index_field()

        from lbmpy.advanced_streaming import Timestep
        from lbmpy.advanced_streaming.indexing import BetweenTimestepsIndexing

        prev_timestep: Timestep = Timestep.BOTH
        streaming_pattern = "pull"

        indexing = BetweenTimestepsIndexing(
            self._pdf_field,
            self._lb_method.stencil,
            prev_timestep,
            streaming_pattern,
            BoundaryIndexType,
            BoundaryIndexType,
        )

        f_out, f_in = indexing.proxy_fields
        dir_symbol = indexing.dir_symbol
        inv_dir = indexing.inverse_dir_symbol

        boundary_assignments = self._bc(
            f_out,
            f_in,
            dir_symbol,
            inv_dir,
            self._lb_method,
            index_field,
            None,  # TODO: Fix force vector
        )
        boundary_assignments = indexing.substitute_proxies(boundary_assignments)

        elements: list[Assignment] = []

        index_arrs_node = indexing.create_code_node()
        elements += index_arrs_node.get_array_declarations()

        if isinstance(self._bc, LbmpyBoundaryWrapper):
            for node in self._bc.lbmpy_object.get_additional_code_nodes(
                self._lb_method
            )[::-1]:
                elements += node.get_array_declarations()

        elements += [Assignment(dir_symbol, index_field[0]("dir"))]
        elements += boundary_assignments.all_assignments

        return elements


class LbmpyBoundaryWrapper(SparseBoundaryDefinition):

    def __init__(self, boundary: LbBoundary):
        self._boundary = boundary
        super().__init__(self._boundary.name)

        self.struct_data = [
            ("x", BoundaryIndexType),
            ("y", BoundaryIndexType),
            ("z", BoundaryIndexType),
            ("dir", BoundaryIndexType),
        ]
        self._idx_struct_type = PsStructType(
            self.struct_data + self._boundary.additional_data,
            self._boundary.name + "Link",
        )
        self._generate_add_data = bool(self._boundary.additional_data)

    def generate_add_data(self):
        return self._generate_add_data

    def index_struct_type(self) -> PsStructType:
        return self._idx_struct_type

    def data_struct_name(self):
        return f"{self._boundary.name}Data"

    def factory_crtp(self) -> str:
        return "walberla::sweepgen::GenericBoundaryFactory"

    def includes(self) -> set[HeaderFile]:
        return {HeaderFile.parse("walberla/sweepgen/boundaries/GenericBoundary.hpp")}

    @property
    def lbmpy_object(self) -> LbBoundary:
        return self._boundary

    def render_link_struct(self, sfg: SfgComposer):
        ctor = sfg.constructor()

        for d_name, d_type in self.struct_data:
            var = SfgVar(d_name, d_type)
            var_ = SfgVar(d_name + "_", d_type)
            ctor.add_param(var_)
            ctor.init(var)(var_)

        if self.generate_add_data():
            data_prop = SfgVar("addData_", PsCustomType(self.data_struct_name()))
            ctor.add_param(data_prop)

            for d_name, d_type in self._boundary.additional_data:
                var = SfgVar(d_name, d_type)
                ctor.init(var)(sfg.expr("{}.{}", data_prop, d_name))

        sfg.struct(self._idx_struct_type.name)(
            *(
                SfgVar(d_name, d_type)
                for d_name, d_type in self.struct_data + self._boundary.additional_data
            ),
            ctor,
        )

    def render_data_struct(self, sfg: SfgComposer):
        if self.generate_add_data():
            ctor = sfg.constructor()
            for d_name, d_type in self._boundary.additional_data:
                var = SfgVar(d_name, d_type)
                var_ = SfgVar(d_name + "_", d_type)
                ctor.add_param(var_)
                ctor.init(var)(var_)

            sfg.struct(self.data_struct_name())(
                *(
                    SfgVar(d_name, d_type)
                    for d_name, d_type in self._boundary.additional_data
                ),
                ctor,
            )

    def __call__(
        self, f_out, f_in, dir_symbol, inv_dir, lb_method, index_field, force_vector
    ) -> Assignment | list[Assignment] | AssignmentCollection:
        return self._boundary(
            f_out, f_in, dir_symbol, inv_dir, lb_method, index_field, force_vector
        )


class FreeSlip(SparseBoundaryDefinition):

    def __init__(self, name: str = "FreeSlip"):
        super().__init__(name)

        self._idx_struct_type = PsStructType(
            (
                ("x", BoundaryIndexType),
                ("y", BoundaryIndexType),
                ("z", BoundaryIndexType),
                ("dir", BoundaryIndexType),
                ("source_offset_x", BoundaryIndexType),
                ("source_offset_y", BoundaryIndexType),
                ("source_offset_z", BoundaryIndexType),
                ("source_dir", BoundaryIndexType),
            ),
            "walberla::sweepgen::FreeSlipLinkInfo",
        )

    def index_struct_type(self):
        return self._idx_struct_type

    def factory_crtp(self) -> str:
        return "walberla::sweepgen::FreeSlipFactory"

    def includes(self) -> set[HeaderFile]:
        return {HeaderFile.parse("walberla/sweepgen/boundaries/FreeSlip.hpp")}

    def __call__(
        self, f_out, f_in, dir_symbol, inv_dir, lb_method, index_field, force_vector
    ) -> Assignment:
        source_cell = (
            index_field("source_offset_x"),
            index_field("source_offset_y"),
            index_field("source_offset_z"),
        )
        source_dir = index_field("source_dir")

        return Assignment(f_in(inv_dir[dir_symbol]), f_out[source_cell](source_dir))
