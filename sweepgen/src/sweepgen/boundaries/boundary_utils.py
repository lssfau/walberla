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

from functools import cache
from abc import abstractmethod

from pystencils import Field, FieldType, TypedSymbol, Assignment
from pystencils.types import PsStructType, create_type

from lbmpy.methods import AbstractLbMethod
from lbmpy.boundaries.boundaryconditions import LbBoundary
from lbmpy.advanced_streaming import Timestep


BoundaryIndexType = create_type("int32")

HbbLinkType = PsStructType(
    [
        ("x", BoundaryIndexType),
        ("y", BoundaryIndexType),
        ("z", BoundaryIndexType),
        ("dir", BoundaryIndexType),
    ],
    "walberla::sweepgen::HbbLink",
)


class WalberlaLbmBoundary:
    idx_struct_type: PsStructType

    @classmethod
    @cache
    def get_index_field(cls):
        return Field(
            "indexVector",
            FieldType.INDEXED,
            cls.idx_struct_type,
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


class GenericHbbWrapper(WalberlaLbmBoundary):
    idx_struct_type = HbbLinkType

    def __init__(self, hbb: LbBoundary):
        self._hbb = hbb

    def boundary_obj(self) -> LbBoundary:
        return self._hbb
