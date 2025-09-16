import pystencils as ps
from lbmpy import Stencil, LBStencil
from pystencilssfg import SourceFileGenerator
from sweepgen.communication import GpuPdfFieldPackInfo
from sweepgen.build_config import DEBUG

DEBUG.use_hip_default()

with SourceFileGenerator() as sfg:
    stencil = LBStencil(Stencil.D3Q19)
    field = ps.fields(f"f({stencil.Q}): double[{stencil.D}D]")
    sfg.generate(GpuPdfFieldPackInfo("PackInfo", stencil, field))
