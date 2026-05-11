from __future__ import annotations

from pystencils import Field, Assignment, AssignmentCollection

from lbmpy.methods import AbstractLbMethod

from pystencilssfg import SfgComposer
from pystencilssfg.composer.custom import CustomGenerator

from ..sweep import Sweep


class GridAlignedNoSlip(CustomGenerator):
    """Zero-velocity half-way bounce back boundary condition at a grid-aligned wall.

    Generate a sweep that applies the zero-velocity half-way bounce back (NoSlip) boundary
    condition to fluid cells touching a wall oriented to the given normal direction.
    The generated sweep must be invoked on the layer of *fluid cells* touching the wall.

    Args:
        name: Name of the generated sweep class
        lb_method: lbmpy Method description
        pdf_field: Symbolic LBM PDF field
        wall_orientation: Vector :math:`\\vec{n} \\in \\{ -1, 0, 1 \\}^d` pointing from the fluid to the wall
    """

    def __init__(
        self,
        name: str,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
        wall_orientation: tuple[int, int, int],
    ):
        self._name = name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._wall_normal = wall_orientation

    def create_sweep(self):
        bc_asm = self._grid_aligned_assignments()
        bc_sweep = Sweep(self._name, bc_asm)
        return bc_sweep

    def generate(self, sfg: SfgComposer):
        sfg.generate(self.create_sweep())

    def _grid_aligned_assignments(self):
        stencil = self._lb_method.stencil

        assert isinstance(self._wall_normal, tuple)
        wall_normal = self._wall_normal
        x_b = (0,) * stencil.D

        f = self._pdf_field

        def wall_aligned(vec):
            """Test if v1 is aligned with v2, i.e. v1 has the same nonzero entries as v2"""
            for x_v, x_n in zip(vec, wall_normal):
                if x_n != 0 and x_v != x_n:
                    return False
            return True

        asms = []

        for i_out, x_w in enumerate(stencil):
            if wall_aligned(x_w):
                i_inverse = stencil.inverse_index(x_w)
                asms.append(Assignment(f[x_w](i_inverse), f[x_b](i_out)))

        return AssignmentCollection(asms)


class GridAlignedFreeSlip(CustomGenerator):
    """Specular-reflection half-way bounce back boundary condition at a grid-aligned wall.

    Generate a sweep that applies a specular-reflection boundary
    condition to fluid cells touching a wall oriented to the given normal direction.
    It can be used to mode both free slip and symmetry boundaries.
    The generated sweep must be invoked on the layer of *fluid cells* touching the wall.

    Args:
        name: Name of the generated sweep class
        lb_method: lbmpy Method description
        pdf_field: Symbolic LBM PDF field
        wall_orientation: Vector :math:`\\vec{n} \\in \\{ -1, 0, 1 \\}^d` pointing from the fluid to the wall
    """

    def __init__(
        self,
        name: str,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
        wall_orientation: tuple[int, int, int],
    ):
        self._name = name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._wall_normal = wall_orientation

    def create_sweep(self):
        bc_asm = self._grid_aligned_assignments()
        bc_sweep = Sweep(self._name, bc_asm)
        return bc_sweep

    def generate(self, sfg: SfgComposer):
        sfg.generate(self.create_sweep())

    def _grid_aligned_assignments(self):
        stencil = self._lb_method.stencil

        assert isinstance(self._wall_normal, tuple)
        wall_normal = self._wall_normal
        inverse_wall_normal = tuple(-c for c in wall_normal)

        f = self._pdf_field

        def is_missing(vec):
            for x_v, x_n in zip(vec, inverse_wall_normal):
                if x_n != 0 and x_v != x_n:
                    return False
            return True

        def wall_mirror(vec):
            return tuple(
                (-x_v if x_n != 0 else x_v) for x_v, x_n in zip(vec, wall_normal)
            )

        asms = []

        for i_missing, c_missing in enumerate(stencil):
            #   Iterate all inbound populations to the fluid cell crossing the wall
            if is_missing(c_missing):
                x_w = tuple(-c for c in c_missing)
                c_refl = wall_mirror(c_missing)
                i_refl = stencil.index(c_refl)
                x_orig = tuple(w - n for w, n in zip(x_w, wall_normal))

                asms.append(Assignment(f[x_w](i_missing), f[x_orig](i_refl)))

        return AssignmentCollection(asms)
