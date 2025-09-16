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
from lbmpy.boundaries import NoSlip, UBB, FixedDensity, SimpleExtrapolationOutflow
from lbmpy import relaxation_rate_from_lattice_viscosity

import sympy as sp
import pystencils as ps

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep, get_build_config
from sweepgen.boundaries import GenericHBB
from sweepgen.symbolic import cell, domain
from sweepgen.prefabs import LbmBulk

from sweepgen.build_config import DEBUG

DEBUG.use_cpu_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    sfg.namespace("FlowAroundSphereExample::gen")
    parser = ArgumentParser()
    parser.add_argument("-t", "--target", choices=["cpu", "gpu"], default="cpu")
    args = parser.parse_args(sfg.context.argv)

    build_cfg = get_build_config(sfg)
    build_cfg.target = ps.Target[args.target.upper()]

    stencil = LBStencil(Stencil.D3Q19)
    omega = sp.symbols("omega")
    inflow_vel = sp.symbols("inflow_vel")

    base_lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.CUMULANT,
        compressible=True,
        relaxation_rate=omega,
    )

    lbm_bulk = LbmBulk(sfg, "LBM", base_lbm_config)
    sfg.generate(lbm_bulk)

    rho, u = lbm_bulk.rho, lbm_bulk.u

    initial_state_assignments = [
        ps.Assignment(rho(), 1),
        ps.Assignment(u(0), inflow_vel)
    ]

    init_fields = Sweep("IntializeMacroFields", initial_state_assignments)
    sfg.generate(init_fields)

    noSlip = GenericHBB(NoSlip(name="NoSlip"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(noSlip)

    inflow_velocity = (inflow_vel, 0, 0)
    ubb = GenericHBB(UBB(inflow_velocity, name="UBB"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(ubb)

    outflow = GenericHBB(SimpleExtrapolationOutflow((1, 0, 0), stencil, name="Outflow"), lbm_bulk.lb_method, lbm_bulk.pdfs)
    sfg.generate(outflow)


