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

from pystencils import (
    Field,
    Assignment,
    AssignmentCollection,
)
from pystencils.types import PsStructType

from lbmpy.methods import AbstractLbMethod
from lbmpy.boundaries.boundaryconditions import LbBoundary

from pystencilssfg import SfgComposer
from pystencilssfg.composer.custom import CustomGenerator
from .boundary_utils import BoundaryIndexType, GenericBoundaryWrapper
from .boundary_utils import WalberlaLbmBoundary
from ..sweep import Sweep
from ..api import SparseIndexList, MemTags
from ..build_config import get_build_config
from ..core.properties import PropertiesContainerBuilder

__all__ = ["FreeSlip"]


class _IrregularSentinel:
    def __repr__(self) -> str:
        return "IRREGULAR"


class GenericLinkwiseBoundary(CustomGenerator):

    IRREGULAR = _IrregularSentinel()


class GenericBoundary(GenericLinkwiseBoundary):
    def __init__(
        self,
        bc: LbBoundary,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
    ):
        self._bc = GenericBoundaryWrapper(bc)
        self._name = bc.name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._stencil = lb_method.stencil

    def generate(self, sfg: SfgComposer) -> None:
        return self._generate_irregular(sfg)

    def _generate_irregular(self, sfg: SfgComposer):
        sfg.include("walberla/sweepgen/boundaries/GenericBoundary.hpp")

        #   Get assignments for bc
        bc_obj = self._bc
        bc_asm = bc_obj.get_assignments(self._lb_method, self._pdf_field)

        #   Build generator config
        bc_cfg = get_build_config(sfg).get_pystencils_config()
        bc_cfg.index_dtype = BoundaryIndexType
        index_field = bc_obj.get_index_field()
        bc_cfg.index_field = index_field

        self._bc.render_data_struct(sfg)
        self._bc.render_link_struct(sfg)

        #   Prepare sweep
        bc_sweep = Sweep(self._name, bc_asm, bc_cfg)

        #   Emit code
        sfg.generate(bc_sweep)

        #   Build factory
        factory_name = f"{self._name}Factory"
        factory_crtp_base = f"walberla::sweepgen::GenericBoundaryFactory< {factory_name} >"
        memtag_t = MemTags.unified if bc_cfg.get_target().is_gpu() else MemTags.host

        index_vector_name = "indexVector"
        index_vector = SparseIndexList(
            bc_obj.idx_struct_type(), memtag_t=memtag_t, ref=True
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
                f"using IdxStruct = {bc_obj.index_struct_name()};",
                f"using DataStruct = {bc_obj.data_struct_name() if bc_obj.generate_add_data() else 'void'};",
                factory_ctor
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

    def _generate_regular(self, sfg: SfgComposer):
        assert False


class NoSlip(GenericLinkwiseBoundary):
    """Zero-velocity half-way bounce back boundary condition.

    Args:
        name: Name of the generated sweep class
        lb_method: lbmpy Method description
        pdf_field: Symbolic LBM PDF field
        wall_orientation: Vector :math:`\\vec{n} \\in \\{ -1, 0, 1 \\}^d` pointing from the fluid to the wall,
            or `NoSlip.IRREGULAR` to indicate a non-grid-aligned boundary
    """

    def __init__(
            self,
            name: str,
            lb_method: AbstractLbMethod,
            pdf_field: Field,
            wall_orientation: tuple[int, int, int] | _IrregularSentinel,
    ):
        self._name = name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._wall_normal = wall_orientation

    def generate(self, sfg: SfgComposer) -> None:
        if self._wall_normal == self.IRREGULAR:
            self._generate_irregular(sfg)
        else:
            self._generate_regular(sfg)

    def _generate_irregular(self, sfg: SfgComposer):
        raise NotImplementedError()

    def _generate_regular(self, sfg: SfgComposer):
        bc_asm = self._grid_aligned_assignments()
        bc_sweep = Sweep(self._name, bc_asm)
        sfg.generate(bc_sweep)

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


class FreeSlip(GenericLinkwiseBoundary):
    """Specular-reflection half-way bounce back boundary condition.

    This boundary handler implements specular reflection of PDFs at a wall located
    at the lattice links' half-way point.
    It can be used to mode both free slip and symmetry boundaries.

    Args:
        name: Name of the generated sweep class
        lb_method: lbmpy Method description
        pdf_field: Symbolic LBM PDF field
        wall_orientation: Vector :math:`\\vec{n} \\in \\{ -1, 0, 1 \\}^d` pointing from the wall into the fluid,
            or `FreeSlip.IRREGULAR` to indicate a non-grid-aligned boundary
    """

    def __init__(
        self,
        name: str,
        lb_method: AbstractLbMethod,
        pdf_field: Field,
        wall_orientation: tuple[int, int, int] | _IrregularSentinel,
    ):
        self._name = name
        self._lb_method = lb_method
        self._pdf_field = pdf_field
        self._wall_normal = wall_orientation

    def generate(self, sfg: SfgComposer) -> None:
        if self._wall_normal == self.IRREGULAR:
            self._generate_irregular(sfg)
        else:
            self._generate_regular(sfg)

    def _generate_irregular(self, sfg: SfgComposer):
        sfg.include("walberla/sweepgen/boundaries/IrregularFreeSlip.hpp")

        #   Get assignments for bc
        bc_obj = WalberlaIrregularFreeSlip()
        bc_asm = bc_obj.get_assignments(self._lb_method, self._pdf_field)

        #   Build generator config
        bc_cfg = get_build_config(sfg).get_pystencils_config()
        bc_cfg.index_dtype = BoundaryIndexType
        index_field = bc_obj.get_index_field()
        bc_cfg.index_field = index_field

        #   Prepare sweep
        bc_sweep = Sweep(self._name, bc_asm, bc_cfg)

        #   Emit code
        sfg.generate(bc_sweep)

        #   Build factory
        factory_name = f"{self._name}Factory"
        factory_crtp_base = (
            f"walberla::sweepgen::IrregularFreeSlipFactory< {factory_name} >"
        )
        memtag_t = MemTags.unified if bc_cfg.get_target().is_gpu() else MemTags.host
        index_vector = SparseIndexList(
            bc_obj.idx_struct_type(), memtag_t=memtag_t, ref=True
        ).var("indexVector")

        sweep_type = bc_sweep.generated_class()
        sweep_ctor_args = {
            f"{self._pdf_field.name}Id": "this->pdfFieldID_",
            f"{index_field.name}": index_vector,
        }

        stencil_name = self._lb_method.stencil.name
        sfg.include(f"stencil/{stencil_name}.h")

        sfg.klass(factory_name, bases=[f"public {factory_crtp_base}"])(
            sfg.public(
                f"using Base = {factory_crtp_base};",
                f"friend class {factory_crtp_base};",
                f"using Stencil = walberla::stencil::{stencil_name};",
                f"using Sweep = {sweep_type.get_dtype().c_string()};",
                f"using memtag_t = {memtag_t.c_string()};" "using Base::Base;",
            ),
            sfg.private(
                sfg.method(
                    "irregularFromIndexVector",
                )
                .returns(sweep_type.get_dtype())
                .inline()(sfg.expr("return {};", sweep_type.ctor(**sweep_ctor_args))),
            ),
        )

    def _generate_regular(self, sfg: SfgComposer):
        bc_asm = self._grid_aligned_assignments()
        bc_sweep = Sweep(self._name, bc_asm)
        sfg.generate(bc_sweep)

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


class WalberlaIrregularFreeSlip(LbBoundary, WalberlaLbmBoundary):

    def __init__(self):
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
            "walberla::sweepgen::IrregularFreeSlipLinkInfo",
        )

    def idx_struct_type(self):
        return self._idx_struct_type

    def boundary_obj(self) -> LbBoundary:
        return self

    def __call__(
        self, f_out, f_in, dir_symbol, inv_dir, lb_method, index_field, force_vector
    ):
        source_cell = (
            index_field("source_offset_x"),
            index_field("source_offset_y"),
            index_field("source_offset_z"),
        )
        source_dir = index_field("source_dir")

        return Assignment(f_in(inv_dir[dir_symbol]), f_out[source_cell](source_dir))
