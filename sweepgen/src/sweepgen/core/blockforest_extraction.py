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

from typing import Callable
import sympy as sp

from pystencils import (
    Assignment,
    AssignmentCollection,
    DEFAULTS
)

from pystencilssfg import SfgComposer
from pystencilssfg.lang import (
    SfgVar,
    AugExpr,
)
from ..api import (
    StructuredBlockForest,
    IBlockPtr,
    CellInterval,
)
from ..symbolic import domain, domain_cell_bb, block, block_cell_bb, cell


class BlockforestParamExtraction:
    def __init__(
        self,
        blockforest: StructuredBlockForest,
        block: IBlockPtr,
    ):
        self._blockforest = blockforest
        self._needs_blockforest = False
        self._block = block
        self._extractions: dict[SfgVar, Callable[[CellInterval | None], AugExpr]] = (
            dict()
        )

    @staticmethod
    def process(asms: AssignmentCollection) -> AssignmentCollection:
        from sympy.core.function import AppliedUndef

        expandable_appls: set[AppliedUndef] = filter(  # type: ignore
            lambda expr: hasattr(expr, "expansion_func"), asms.atoms(AppliedUndef)
        )

        subs: dict[AppliedUndef, sp.Symbol] = dict()
        for appl in expandable_appls:
            expansion: sp.Expr = appl.expansion_func(*appl.args)  # type: ignore
            symb = next(asms.subexpression_symbol_generator)  # type: ignore
            asms.subexpressions.insert(0, Assignment(symb, expansion))
            subs[appl] = symb

        return asms.new_with_substitutions(subs)

    def filter_params(self, params: set[SfgVar]) -> tuple[set[SfgVar], bool]:
        params_filtered = set()

        def get_mapping(param: SfgVar, coord: int):
            #   The actual extractions are modelled as lambdas
            #   since the cell interval is not known at this point
            #   Also, we differentiate between extractions
            #   that require the blockforest object and those which don't,
            #   such that the blockforest pointer only needs to be added to the sweep
            #   class if it is really required

            without_blockforest = {
                block.aabb_min[coord].name: lambda _: (
                    self._block.getAABB().min()[coord]
                ),
                block.aabb_max[coord].name: lambda _: (
                    self._block.getAABB().max()[coord]
                ),
                block.ci_min[coord].name: lambda ci: (
                    AugExpr.format("0") if ci is None else ci.min()[coord]
                ),
                block.ci_max[coord].name: lambda ci: (
                    AugExpr.format("0") if ci is None else ci.max()[coord]
                ),
                #   If spatial counters escape the kernel (e.g. the z-counter in a 2D-kernel),
                #   they are set to zero
                DEFAULTS.spatial_counters[coord].name: lambda _: AugExpr.format("0")
            }

            with_blockforest = {
                domain.aabb_min[coord].name: lambda _: (
                    self._blockforest.getDomain().min()[coord]
                ),
                domain.aabb_max[coord].name: lambda _: (
                    self._blockforest.getDomain().max()[coord]
                ),
                domain_cell_bb.cell_bb_min[coord].name: lambda _: (
                    self._blockforest.getDomainCellBB(
                        self._blockforest.getLevel(self._block.deref())
                    ).min()[coord]
                ),
                domain_cell_bb.cell_bb_max[coord].name: lambda _: (
                    self._blockforest.getDomainCellBB(
                        self._blockforest.getLevel(self._block.deref())
                    ).max()[coord]
                ),
                block_cell_bb.cell_bb_min[coord].name: lambda _: (
                    self._blockforest.getBlockCellBB(self._block.deref()).min()[coord]
                ),
                block_cell_bb.cell_bb_max[coord].name: lambda _: (
                    self._blockforest.getBlockCellBB(self._block.deref()).max()[coord]
                ),
                cell.cell_extents[coord].name: lambda _: (
                    self._blockforest.cell_extents(
                        coord, self._blockforest.getLevel(self._block.deref())
                    )
                ),
            }

            if param.name in without_blockforest:
                return without_blockforest[param.name], False
            elif param.name in with_blockforest:
                return with_blockforest[param.name], True
            else:
                return None, False

        needs_blockforest = False
        for param in params:
            for coord in range(3):
                extraction, bfs = get_mapping(param, coord)
                if extraction is not None:
                    self._extractions[param] = extraction
                    needs_blockforest = needs_blockforest or bfs
                    break
            else:
                params_filtered.add(param)

        return params_filtered, needs_blockforest

    def render_extractions(self, sfg: SfgComposer, ci: CellInterval | None):
        return tuple(sfg.init(p)(expr(ci)) for p, expr in self._extractions.items())
