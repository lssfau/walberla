import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_nonuniform_pack_info_for_field

with CodeGeneration() as ctx:
    layout = 'fzyx'
    field = ps.fields("field(1): float64[3D]", layout=layout)
    generate_nonuniform_pack_info_for_field(ctx, 'ScalarFieldNonUniformCommunicationCPU', field, target=ps.Target.CPU)