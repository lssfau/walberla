import operator as op
import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_pack_info_for_field

with CodeGeneration() as ctx:

    layout = 'fzyx'

    field = ps.fields("field: int32[3D]", layout=layout)

    # communication
    generate_pack_info_for_field(ctx, 'ScalarFieldCommunicationGPU', field, target='gpu')
    generate_pack_info_for_field(ctx, 'ScalarFieldPullReductionGPU', field, target='gpu', operator=op.add,
                                 gl_to_inner=True)
