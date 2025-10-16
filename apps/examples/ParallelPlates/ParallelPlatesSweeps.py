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

from argparse import ArgumentParser
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

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.boundaries import GenericBoundary
from sweepgen.symbolic import cell, domain
from sweepgen.prefabs import LbmBulk

from sweepgen.build_config import DEBUG

DEBUG.use_cpu_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    sfg.namespace("ParallelPlates::gen")
    parser = ArgumentParser()
    parser.add_argument("-t", "--target", choices=["cpu", "gpu"], default="cpu")
    args = parser.parse_args(sfg.context.argv)

    build_cfg = get_build_config(sfg)
    build_cfg.target = ps.Target[args.target.upper()]

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
        r = sp.Symbol("r")

        poiseuille_analytical_asms = [
            ps.Assignment(r, sp.Abs(cell.z() - R)),
            ps.Assignment(rho(), 1),
            ps.Assignment(u(0), a_x / (2 * nu) * (R**2 - r**2)),
            ps.Assignment(u(1), 0),
            ps.Assignment(u(2), 0),
        ]
        sfg.generate(Sweep("SetAnalyticalSolution", poiseuille_analytical_asms))

        ux = sp.Symbol("ux")
        error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

        poiseuille_velocity_error_lmax = [
            ps.Assignment(r, sp.Abs(cell.z() - R)),
            ps.Assignment(ux, a_x / (2 * nu) * (R**2 - r**2)),
            ps.MaxReductionAssignment(error_ux, sp.Abs(u(0) - ux))
        ]
        sfg.generate(Sweep("VelocityErrorLmax", poiseuille_velocity_error_lmax))

    #   Setup for Couette Flow
    with sfg.namespace("Couette"):
        lbm_bulk = LbmBulk(sfg, "LBM", base_lbm_config)
        sfg.generate(lbm_bulk)

        rho, u = lbm_bulk.rho, lbm_bulk.u

        couette_analytical_asms = [
            ps.Assignment(rho(), 1),
            ps.Assignment(u(0), u_max * cell.z() / domain.z_max()),
            ps.Assignment(u(1), 0),
            ps.Assignment(u(2), 0),
        ]
        sfg.generate(Sweep("SetAnalyticalSolution", couette_analytical_asms))

        ux = sp.Symbol("ux")
        error_ux = ps.TypedSymbol("error_ux", ps.DynamicType.NUMERIC_TYPE)

        couette_velocity_error_lmax = [
            ps.Assignment(ux, u_max * cell.z() / domain.z_max()),
            ps.MaxReductionAssignment(error_ux, sp.Abs(u(0) - ux))
        ]
        sfg.generate(Sweep("VelocityErrorLmax", couette_velocity_error_lmax))

    #   Boundary Conditions

    noSlip = GenericBoundary(NoSlip(name="NoSlip"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(noSlip)

    wall_velocity = (u_max, 0, 0)
    ubb = GenericBoundary(UBB(wall_velocity, name="UBB"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(ubb)
