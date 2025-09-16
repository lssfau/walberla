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

import sympy as sp

import inspect

from pystencils import TypedSymbol, DynamicType, tcast
from pystencils.defaults import DEFAULTS


def expandable(name: str):
    def wrap(expansion_func):
        nargs = len(inspect.signature(expansion_func).parameters)
        return sp.Function(
            name, nargs=nargs, expansion_func=staticmethod(expansion_func)
        )

    return wrap


class DomainCoordinates:
    aabb_min = tuple(
        TypedSymbol(f"blockforest_aabb_min_{i}", DynamicType.NUMERIC_TYPE)
        for i in range(3)
    )

    aabb_max = tuple(
        TypedSymbol(f"blockforest_aabb_max_{i}", DynamicType.NUMERIC_TYPE)
        for i in range(3)
    )

    @expandable("domain.x_min")
    @staticmethod
    def x_min():
        """The minimum X coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_min[0]

    @expandable("domain.y_min")
    @staticmethod
    def y_min():
        """The minimum Y coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_min[1]

    @expandable("domain.z_min")
    @staticmethod
    def z_min():
        """The minimum Z coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_min[2]

    @expandable("domain.x_max")
    @staticmethod
    def x_max():
        """The maximum X coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_max[0]

    @expandable("domain.y_max")
    @staticmethod
    def y_max():
        """The maximum Y coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_max[1]

    @expandable("domain.z_max")
    @staticmethod
    def z_max():
        """The maximum Z coordinate of the domain axis-aligned bounding box"""
        return DomainCoordinates.aabb_max[2]


domain = DomainCoordinates


class DomainCellBoundingBox:
    """Represents the domain cell bounding box on the current block's refinement level.

    See ``StructuredBlockStorage::getDomainCellBB``.
    """

    cell_bb_min = tuple(
        TypedSymbol(f"domain_cell_b_min_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )

    cell_bb_max = tuple(
        TypedSymbol(f"domain_cell_b_max_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )

    @expandable("domain_cell_bb.x_min")
    @staticmethod
    def x_min():
        """The X coordinate of the domain's minimum cell index"""
        return DomainCellBoundingBox.cell_bb_min[0]

    @expandable("domain_cell_bb.y_min")
    @staticmethod
    def y_min():
        """The Y coordinate of the domain's minimum cell index"""
        return DomainCellBoundingBox.cell_bb_min[1]

    @expandable("domain_cell_bb.z_min")
    @staticmethod
    def z_min():
        """The Z coordinate of the domain's minimum cell index"""
        return DomainCellBoundingBox.cell_bb_min[2]

    @expandable("domain_cell_bb.x_max")
    @staticmethod
    def x_max():
        """The X coordinate of the domain's maximum cell index"""
        return DomainCellBoundingBox.cell_bb_max[0]

    @expandable("domain_cell_bb.y_max")
    @staticmethod
    def y_max():
        """The Y coordinate of the domain's maximum cell index"""
        return DomainCellBoundingBox.cell_bb_max[1]

    @expandable("domain_cell_bb.z_max")
    @staticmethod
    def z_max():
        """The Z coordinate of the domain's maximum cell index"""
        return DomainCellBoundingBox.cell_bb_max[2]


domain_cell_bb = DomainCellBoundingBox


class BlockCoordinates:
    aabb_min = tuple(
        TypedSymbol(f"block_aabb_min_{i}", DynamicType.NUMERIC_TYPE) for i in range(3)
    )
    aabb_max = tuple(
        TypedSymbol(f"block_aabb_max_{i}", DynamicType.NUMERIC_TYPE) for i in range(3)
    )
    ci_min = tuple(
        TypedSymbol(f"cell_interval_min_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )
    ci_max = tuple(
        TypedSymbol(f"cell_interval_max_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )


block = BlockCoordinates


class BlockCellBoundingBox:
    """Represents the cell bounding box of the current block.

    See ``StructuredBlockStorage::getBlockCellBB``.
    """

    cell_bb_min = tuple(
        TypedSymbol(f"block_cell_bb_min_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )

    cell_bb_max = tuple(
        TypedSymbol(f"block_cell_bb_max_{i}", DynamicType.INDEX_TYPE) for i in range(3)
    )

    @expandable("block_cell_bb.x_min")
    @staticmethod
    def x_min():
        """The X coordinate of the current block's minimum cell index"""
        return BlockCellBoundingBox.cell_bb_min[0]

    @expandable("block_cell_bb.y_min")
    @staticmethod
    def y_min():
        """The Y coordinate of the current block's minimum cell index"""
        return BlockCellBoundingBox.cell_bb_min[1]

    @expandable("block_cell_bb.z_min")
    @staticmethod
    def z_min():
        """The Z coordinate of the current block's minimum cell index"""
        return BlockCellBoundingBox.cell_bb_min[2]

    @expandable("block_cell_bb.x_max")
    @staticmethod
    def x_max():
        """The X coordinate of the current block's maximum cell index"""
        return BlockCellBoundingBox.cell_bb_max[0]

    @expandable("block_cell_bb.y_max")
    @staticmethod
    def y_max():
        """The Y coordinate of the current block's maximum cell index"""
        return BlockCellBoundingBox.cell_bb_max[1]

    @expandable("block_cell_bb.z_max")
    @staticmethod
    def z_max():
        """The Z coordinate of the current block's maximum cell index"""
        return BlockCellBoundingBox.cell_bb_max[2]


block_cell_bb = BlockCellBoundingBox


class CellCoordinates:
    cell_extents = tuple(
        TypedSymbol(f"cell_extents_{i}", DynamicType.NUMERIC_TYPE) for i in range(3)
    )

    @expandable("cell.x")
    @staticmethod
    def x():
        """X coordinate of the current cell's center in the global coordinate system"""
        return BlockCoordinates.aabb_min[0] + CellCoordinates.cell_extents[0] * (
            tcast.as_numeric(
                block.ci_min[0] + DEFAULTS.spatial_counters[0]
            )
            + sp.Rational(1, 2)
        )

    @expandable("cell.y")
    @staticmethod
    def y():
        """Y coordinate of the current cell's center in the global coordinate system"""
        return BlockCoordinates.aabb_min[1] + CellCoordinates.cell_extents[1] * (
            tcast.as_numeric(
                block.ci_min[1] + DEFAULTS.spatial_counters[1]
            )
            + sp.Rational(1, 2)
        )

    @expandable("cell.z")
    @staticmethod
    def z():
        """Z coordinate of the current cell's center in the global coordinate system"""
        return BlockCoordinates.aabb_min[2] + CellCoordinates.cell_extents[2] * (
            tcast.as_numeric(
                block.ci_min[2] + DEFAULTS.spatial_counters[2]
            )
            + sp.Rational(1, 2)
        )

    @expandable("cell.local_x")
    @staticmethod
    def local_x():
        """X coordinate of the current cell's center in the current block's local coordinate system"""
        return CellCoordinates.cell_extents[0] * (
            tcast.as_numeric(block.ci_min[0] + DEFAULTS.spatial_counters[0])
            + sp.Rational(1, 2)
        )

    @expandable("cell.local_y")
    @staticmethod
    def local_y():
        """Y coordinate of the current cell's center in the current block's local coordinate system"""
        return CellCoordinates.cell_extents[1] * (
            tcast.as_numeric(block.ci_min[1] + DEFAULTS.spatial_counters[1])
            + sp.Rational(1, 2)
        )

    @expandable("cell.local_z")
    @staticmethod
    def local_z():
        """Z coordinate of the current cell's center in the current block's local coordinate system"""
        return CellCoordinates.cell_extents[2] * (
            tcast.as_numeric(block.ci_min[2] + DEFAULTS.spatial_counters[2])
            + sp.Rational(1, 2)
        )

    @expandable("cell.dx")
    @staticmethod
    def dx():
        """Size of this cell in x-direction"""
        return CellCoordinates.cell_extents[0]

    @expandable("cell.dy")
    @staticmethod
    def dy():
        """Size of this cell in y-direction"""
        return CellCoordinates.cell_extents[1]

    @expandable("cell.dz")
    @staticmethod
    def dz():
        """Size of this cell in z-direction"""
        return CellCoordinates.cell_extents[2]


cell = CellCoordinates


class CellIndex:
    """Represents the index of the current cell in the global and block-local cell grids."""

    @expandable("cell_index.x_global")
    @staticmethod
    def x_global():
        """X component of the current cell's index in the global cell grid on the current block's refinement level."""
        return (
            BlockCellBoundingBox.cell_bb_min[0]
            + block.ci_min[0]
            + DEFAULTS.spatial_counters[0]
        )

    @expandable("cell_index.y_global")
    @staticmethod
    def y_global():
        """Y component of the current cell's index in the global cell grid on the current block's refinement level."""
        return (
            BlockCellBoundingBox.cell_bb_min[1]
            + block.ci_min[1]
            + DEFAULTS.spatial_counters[1]
        )

    @expandable("cell_index.z_global")
    @staticmethod
    def z_global():
        """Z component of the current cell's index in the global cell grid on the current block's refinement level."""
        return (
            BlockCellBoundingBox.cell_bb_min[2]
            + block.ci_min[2]
            + DEFAULTS.spatial_counters[2]
        )

    @expandable("cell_index.x_local")
    @staticmethod
    def x_local():
        """X component of the current cell's index in the current block's local cell grid."""
        return block.ci_min[0] + DEFAULTS.spatial_counters[0]

    @expandable("cell_index.y_local")
    @staticmethod
    def y_local():
        """Y component of the current cell's index in the current block's local cell grid."""
        return block.ci_min[1] + DEFAULTS.spatial_counters[1]

    @expandable("cell_index.z_local")
    @staticmethod
    def z_local():
        """Z component of the current cell's index in the current block's local cell grid."""
        return block.ci_min[2] + DEFAULTS.spatial_counters[2]


cell_index = CellIndex
