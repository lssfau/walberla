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

from pystencils import fields, CreateKernelConfig, Field

from pystencilssfg import SfgComposer, SfgContext
from pystencilssfg.composer import SfgIComposer
from pystencilssfg.composer.custom import CustomGenerator
from lbmpy import LBMConfig, LBMOptimisation, create_lb_update_rule
from lbmpy.methods import LbmCollisionRule, AbstractLbMethod
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

from ..build_config import get_build_config
from ..sweep import Sweep


class LbmBulk(CustomGenerator):
    """Generate a sweep group for lattice Boltzmann bulk dynamics.

    Generate a sweep group for a given lattice Boltzmann method definition,
    including:

    - A PDF initialization sweep
    - A stream/collide sweep

    **Example:**

    .. code-block:: python

        lb_config = LBMConfig(
            stencil=Stencil.D3Q19,
            method=Method.CENTRAL_MOMENT,
            relaxation_rate=sp.Symbol("omega")
        )
        lb_bulk = LbmBulk(sfg, "MyLBM", lb_config)
        sfg.generate(lb_bulk)

    .. admonition:: Example Apps
        :class: seealso

        - :walberla-example-app:`DoubleShearLayer`
        - :walberla-example-app:`ParallelPlates`

    Args:
        sfg: An SFG composer or context that holds a waLBerla build configuration
        name: Name of the sweep group
        lbm_config: lbmpy method configuration defining the lattice Boltzmann method
        lbm_optimisation: Optimization settings for lbmpy
        gen_config: pystencils code generator configuration
    """

    def __init__(
        self,
        sfg: SfgIComposer | SfgContext,
        name: str,
        lbm_config: LBMConfig | None = None,
        lbm_optimisation: LBMOptimisation | None = None,
        gen_config: CreateKernelConfig | None = None,
    ):
        self._name = name
        if lbm_config is None:
            lbm_config = LBMConfig()

        if lbm_config.streaming_pattern != "pull":
            raise ValueError("Only the `pull` streaming pattern is currently supported.")

        self._lbm_config = lbm_config

        base_gen_config = get_build_config(sfg).get_pystencils_config()
        if gen_config is not None:
            base_gen_config.override(gen_config)
        self._gen_config = base_gen_config

        self._stencil = lbm_config.stencil
        Q, D = self._stencil.Q, self._stencil.D

        self._pdfs: Field
        self._pdfs_tmp: Field
        self._rho: Field
        self._u: Field

        dtype = self._gen_config.get_option("default_dtype")

        if lbm_optimisation is None:
            lbm_optimisation = LBMOptimisation(
                simplification="auto",
                cse_global=False,
                field_layout="fzyx",
            )

        if lbm_optimisation.symbolic_field is None:
            self._pdfs = lbm_optimisation.symbolic_field = Field.create_generic(
                "f",
                spatial_dimensions=D,
                index_shape=(Q,),
                layout=lbm_optimisation.field_layout,
                dtype=dtype,
            )

        if lbm_optimisation.symbolic_temporary_field is None:
            self._pdfs_tmp = lbm_optimisation.symbolic_temporary_field = (
                Field.create_generic(
                    "f_tmp",
                    spatial_dimensions=D,
                    index_shape=(Q,),
                    layout=lbm_optimisation.field_layout,
                    dtype=dtype,
                )
            )

        self._rho, self._u = fields(  # type: ignore
            f"rho, u({D}): {dtype}[{D}D]", layout=lbm_optimisation.field_layout
        )

        output_fields = dict({"density": self._rho, "velocity": self._u})

        self._update_rule: LbmCollisionRule = create_lb_update_rule(
            lbm_config=lbm_config, lbm_optimisation=lbm_optimisation, output=output_fields
        )

        self._init_rule = macroscopic_values_setter(
            lb_method=self._update_rule.method,
            density=self._rho.center,
            velocity=self._u.center_vector,
            pdfs=self._pdfs,
            set_pre_collision_pdfs=True,
        )

    @property
    def pdfs(self) -> Field:
        """The symbolic PDF field"""
        return self._pdfs

    @property
    def rho(self) -> Field:
        """The symbolic density field"""
        return self._rho

    @property
    def u(self) -> Field:
        """The symbolic velocity field"""
        return self._u

    @property
    def update_rule(self) -> LbmCollisionRule:
        """The symbolic LB update rule"""
        return self._update_rule

    @property
    def lb_method(self) -> AbstractLbMethod:
        """The LB method object"""
        return self._update_rule.method

    def generate(self, sfg: SfgComposer) -> None:
        with sfg.namespace(self._name):
            stream_collide_sweep = Sweep(
                "StreamCollide", self._update_rule, self._gen_config
            )
            stream_collide_sweep.swap_fields(self._pdfs, self._pdfs_tmp)
            sfg.generate(stream_collide_sweep)

            init_sweep = Sweep("InitPdfs", self._init_rule, self._gen_config)
            sfg.generate(init_sweep)
