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
import sweepgen as sg
from sweepgen.build_config import DEBUG
from sweepgen.symbolic import cell
from sweepgen.prefabs import LbmBulk

DEBUG.use_cpu_default()


with SourceFileGenerator() as sfg:
    sg.Sweep.use_v8core_fields()
    sfg.namespace("DoubleShearLayer::gen")

    stencil = LBStencil(Stencil.D3Q19)
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,
        compressible=True,
        relaxation_rate=sp.Symbol("omega"),
    )

    lbm_bulk = LbmBulk(sfg, "LBM", lbm_config)
    sfg.generate(lbm_bulk)

    rho, u = lbm_bulk.rho, lbm_bulk.u

    #   Initial State

    @sg.flow.generate_sweep(sfg)
    def SetInitialState(_eq):
        u_0, kappa, delta = sp.symbols("u_0, kappa, delta")

        _eq.store[rho()] = 1

        _eq.store[u(0)] = sp.Piecewise(
            (
                u_0 * sp.tanh(kappa * (cell.y() - sp.Rational(1, 4))),
                cell.y() <= sp.Rational(1, 2),
            ),
            (
                u_0 * sp.tanh(kappa * (sp.Rational(3, 4) - cell.y())),
                cell.y() > sp.Rational(1, 2),
            ),
        )

        _eq.store[u(1)] = (
            delta * u_0 * sp.sin(2 * sp.pi * (cell.x() + sp.Rational(1, 4)))
        )
        _eq.store[u(2)] = 0

    #   end initial state

    #   Compute Vorticity

    @sg.flow.generate_sweep(sfg)
    def ComputeVorticity(_eq):
        dvx, duy = sp.symbols("dvx, duy")
        vorticity = ps.fields(f"vorticity: double[{stencil.D}D]", layout="fzyx")

        _eq.let[dvx] = (u[1, 0, 0](1) - u[-1, 0, 0](1)) / (2 * cell.dx())
        _eq.let[duy] = (u[0, 1, 0](0) - u[0, -1, 0](0)) / (2 * cell.dy())
        _eq.store[vorticity()] = (dvx - duy) / 2

    #   end vorticity
