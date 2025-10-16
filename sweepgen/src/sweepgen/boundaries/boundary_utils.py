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

from abc import abstractmethod

from pystencils import Field, FieldType, TypedSymbol, Assignment
from pystencils.types import PsStructType, create_type, PsCustomType
from pystencilssfg import SfgComposer
from pystencilssfg.lang import SfgVar

from lbmpy.methods import AbstractLbMethod
from lbmpy.boundaries.boundaryconditions import LbBoundary
from lbmpy.advanced_streaming import Timestep

BoundaryIndexType = create_type("int32")


class WalberlaLbmBoundary:

    def get_index_field(self):
        return Field(
            "indexVector",
            FieldType.INDEXED,
            self.idx_struct_type(),
            (0,),
            (TypedSymbol("indexVectorLength", BoundaryIndexType), 1),
            (1, 1),
        )

    @abstractmethod
    def boundary_obj(self) -> LbBoundary: ...

    def get_assignments(
        self,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
        prev_timestep: Timestep = Timestep.BOTH,
        streaming_pattern="pull",
    ) -> list[Assignment]:
        index_field = self.get_index_field()

        from lbmpy.advanced_streaming.indexing import BetweenTimestepsIndexing

        indexing = BetweenTimestepsIndexing(
            pdf_field,
            lb_method.stencil,
            prev_timestep,
            streaming_pattern,
            BoundaryIndexType,
            BoundaryIndexType,
        )

        f_out, f_in = indexing.proxy_fields
        dir_symbol = indexing.dir_symbol
        inv_dir = indexing.inverse_dir_symbol

        boundary_assignments = self.boundary_obj()(
            f_out,
            f_in,
            dir_symbol,
            inv_dir,
            lb_method,
            index_field,
            None,  # TODO: Fix force vector
        )
        boundary_assignments = indexing.substitute_proxies(boundary_assignments)

        elements: list[Assignment] = []

        index_arrs_node = indexing.create_code_node()
        elements += index_arrs_node.get_array_declarations()

        for node in self.boundary_obj().get_additional_code_nodes(lb_method)[::-1]:
            elements += node.get_array_declarations()

        elements += [Assignment(dir_symbol, index_field[0]("dir"))]
        elements += boundary_assignments.all_assignments

        return elements


class GenericBoundaryWrapper(WalberlaLbmBoundary):

    def __init__(self, boundary: LbBoundary):
        self._boundary = boundary

        self.struct_data = [
            ("x", BoundaryIndexType),
            ("y", BoundaryIndexType),
            ("z", BoundaryIndexType),
            ("dir", BoundaryIndexType),
        ]
        self._idx_struct_type = PsStructType(
            self.struct_data + self._boundary.additional_data,
            self._boundary.name+"Link")
        self._generate_add_data = bool(self._boundary.additional_data)

    def generate_add_data(self):
        return self._generate_add_data

    def idx_struct_type(self) -> PsStructType:
        return self._idx_struct_type

    def boundary_obj(self) -> LbBoundary:
        return self._boundary

    def index_struct_name(self):
        return f"{self._boundary.name}Link"

    def data_struct_name(self):
        return f"{self._boundary.name}Data"

    def render_link_struct(self, sfg: SfgComposer):
        ctor = sfg.constructor()

        for d_name, d_type in self.struct_data:
            var = SfgVar(d_name, d_type)
            var_ = SfgVar(d_name+"_", d_type)
            ctor.add_param(var_)
            ctor.init(var)(var_)

        if self.generate_add_data():
            data_prop = SfgVar("addData_", PsCustomType(self.data_struct_name()))
            ctor.add_param(data_prop)

            for d_name, d_type in self._boundary.additional_data:
                var = SfgVar(d_name, d_type)
                ctor.init(var)(sfg.expr("{}.{}", data_prop, d_name))

        sfg.struct(self.index_struct_name())(
            *(SfgVar(d_name, d_type) for d_name, d_type in self.struct_data + self._boundary.additional_data),
            ctor,
        )

    def render_data_struct(self, sfg: SfgComposer):
        if self.generate_add_data():
            ctor = sfg.constructor()
            for d_name, d_type in self._boundary.additional_data:
                var = SfgVar(d_name, d_type)
                var_ = SfgVar(d_name+"_", d_type)
                ctor.add_param(var_)
                ctor.init(var)(var_)

            sfg.struct(self.data_struct_name())(
                *(SfgVar(d_name, d_type) for d_name, d_type in self._boundary.additional_data),
                ctor,
            )
