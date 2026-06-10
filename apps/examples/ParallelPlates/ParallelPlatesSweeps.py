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

from dataclasses import replace

from lbmpy import (
    LBStencil,
    Stencil,
    LBMConfig,
    Method,
    ForceModel,
)
from lbmpy.boundaries import NoSlip, UBB
from lbmpy import relaxation_rate_from_lattice_viscosity

import sympy as sp
import pystencils as ps
import sweepgen as sg

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep
from sweepgen.boundaries import SparseBoundary
from sweepgen.symbolic import cell, domain
from sweepgen.prefabs import LbmBulk

from sweepgen.build_config import DEBUG

DEBUG.use_cuda_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    sfg.namespace("ParallelPlates::gen")
    Sweep.use_v8core_fields()

    stencil = LBStencil(Stencil.D3Q19)
    nu, u_max, rho = sp.symbols("nu, u_max, rho")

    base_lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.TRT,
        compressible=True,
        relaxation_rate=relaxation_rate_from_lattice_viscosity(nu),
    )

    #   Setup for Poiseuille Flow
    with sfg.namespace("Poiseuille"):
        R = domain.z_max() / 2
        a_x = 2 * u_max * nu / R**2

        lbm_config = replace(
            base_lbm_config,
            force_model=ForceModel.GUO,
            force=(a_x * rho, 0, 0),
        )

        lbm_bulk = LbmBulk(sfg, "LBM", lbm_config)
        sfg.generate(lbm_bulk)

        rho, u = lbm_bulk.rho, lbm_bulk.u

        @sg.flow.generate_sweep(sfg)
        def SetAnalyticalSolution(_eq):
            r = sp.Symbol("r")

            _eq.let[r] = sp.Abs(cell.z() - R)

            _eq.store[rho()] = 1
            _eq.store[u(0)] = a_x / (2 * nu) * (R**2 - r**2)
            _eq.store[u(1)] = 0
            _eq.store[u(2)] = 0

        @sg.flow.generate_sweep(sfg)
        def VelocityErrorLmax(_eq):
            r, ux = sp.symbols("r, ux")
            error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

            _eq.let[r] = sp.Abs(cell.z() - R)
            _eq.let[ux] = a_x / (2 * nu) * (R**2 - r**2)
            _eq.reduce[error_ux, "max"] = sp.Abs(u(0) - ux)

    #   Setup for Couette Flow
    with sfg.namespace("Couette"):
        lbm_bulk = LbmBulk(sfg, "LBM", base_lbm_config)
        sfg.generate(lbm_bulk)

        rho, u = lbm_bulk.rho, lbm_bulk.u

        @sg.flow.generate_sweep(sfg)
        def SetAnalyticalSolution(_eq):
            _eq.store[rho()] = 1
            _eq.store[u(0)] = u_max * cell.z() / domain.z_max()
            _eq.store[u(1)] = 0
            _eq.store[u(2)] = 0

        @sg.flow.generate_sweep(sfg)
        def VelocityErrorLmax(_eq):
            ux = sp.Symbol("ux")
            error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

            _eq.let[ux] = u_max * cell.z() / domain.z_max()
            _eq.reduce[error_ux, "max"] = sp.Abs(u(0) - ux)

    #   Boundary Conditions

    noSlip = SparseBoundary(NoSlip(name="NoSlip"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(noSlip)

    wall_velocity = (u_max, 0, 0)
    ubb = SparseBoundary(
        UBB(wall_velocity, name="UBB"), lbm_bulk.lb_method, lbm_bulk.pdfs
    )
    sfg.generate(ubb)
