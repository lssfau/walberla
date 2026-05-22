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

from pystencilssfg import SourceFileGenerator
from sweepgen import Sweep
from sweepgen.prefabs import LbmBulk

with SourceFileGenerator() as sfg:
    Sweep.use_v8core_fields()
    sfg.namespace("TestCommunicationHiding::gen")

    stencil = LBStencil(Stencil.D3Q19)
    lbm_config = LBMConfig(
        stencil=stencil,
        method=Method.SRT,
        compressible=True,
        relaxation_rate=sp.Symbol("omega"),
    )

    lbm_bulk = LbmBulk(sfg, "LBM", lbm_config)
    sfg.generate(lbm_bulk)


