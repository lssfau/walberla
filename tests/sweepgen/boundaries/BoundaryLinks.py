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

from lbmpy import LBStencil, Stencil, LBMConfig, Method, create_lb_method
from lbmpy.boundaries import NoSlip, QuadraticBounceBack

import sympy as sp
import pystencils as ps

from pystencilssfg import SourceFileGenerator
from sweepgen.boundaries import GenericBoundary
from sweepgen.prefabs import LbmBulk
from sweepgen import get_build_config

from sweepgen.build_config import DEBUG
DEBUG.use_cpu_default()

with SourceFileGenerator(keep_unknown_argv=True) as sfg:
    sfg.namespace("BoundaryLinks::gen")
    build_cfg = get_build_config(sfg)
    build_cfg.target = ps.Target.GenericCPU

    stencil = LBStencil(Stencil.D3Q7)
    omega = sp.symbols("omega")

    dtype = build_cfg.get_pystencils_config().get_option("default_dtype")
    pdfs = ps.fields( f"pdfs({stencil.Q}): {dtype}[{stencil.D}D]", layout="fzyx")

    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,
        compressible=True,
        relaxation_rate=omega,
    )

    lb_method = create_lb_method(lbm_config)

    noSlip = GenericBoundary(NoSlip(name="NoSlip"), lb_method, pdfs)
    sfg.generate(noSlip)

    qbb = GenericBoundary(QuadraticBounceBack(omega, name="QBB"), lb_method, pdfs)
    sfg.generate(qbb)


