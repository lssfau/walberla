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

from lbmpy import (
    LBStencil,
    Stencil,
    LBMConfig,
    Method,
)

import sympy as sp
import pystencils as ps

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.symbolic import cell
from sweepgen.prefabs import LbmBulk


with SourceFileGenerator() as sfg:
    sfg.namespace("DoubleShearLayer::gen")

    get_build_config(sfg).target = ps.Target.CurrentCPU

    stencil = LBStencil(Stencil.D3Q19)
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,
        compressible=True,
        relaxation_rate=sp.Symbol("omega"),
    )

    lbm_bulk = LbmBulk(sfg, "LBM", lbm_config)
    sfg.generate(lbm_bulk)

    #   Initial State
    rho, u = lbm_bulk.rho, lbm_bulk.u
    u_0, kappa, delta = sp.symbols("u_0, kappa, delta")

    initial_state_assignments = [
        ps.Assignment(rho(), 1),
        ps.Assignment(
            u(0),
            sp.Piecewise(
                (
                    u_0 * sp.tanh(kappa * (cell.y() - sp.Rational(1, 4))),
                    cell.y() <= sp.Rational(1, 2),
                ),
                (
                    u_0 * sp.tanh(kappa * (sp.Rational(3, 4) - cell.y())),
                    cell.y() > sp.Rational(1, 2),
                ),
            ),
        ),
        ps.Assignment(
            u(1), delta * u_0 * sp.sin(2 * sp.pi * (cell.x() + sp.Rational(1, 4)))
        ),
        ps.Assignment(u(2), 0),
    ]

    initial_state_sweep = Sweep("SetInitialState", initial_state_assignments)
    sfg.generate(initial_state_sweep)

    #   Compute Vorticity

    dvx, duy = sp.symbols("dvx, duy")
    vorticity = ps.fields(f"vorticity: double[{stencil.D}D]", layout="fzyx")
    vorticity_assignments = [
        ps.Assignment(dvx, (u[1, 0, 0](1) - u[-1, 0, 0](1)) / (2 * cell.dx())),
        ps.Assignment(duy, (u[0, 1, 0](0) - u[0, -1, 0](0)) / (2 * cell.dy())),
        ps.Assignment(vorticity(), (dvx - duy) / 2),
    ]
    sfg.generate(Sweep("ComputeVorticity", vorticity_assignments))
