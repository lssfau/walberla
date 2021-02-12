from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_update_rule
from lbmpy.advanced_streaming import Timestep
from lbmpy.stencils import get_stencil
from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla.packinfo import generate_lb_pack_info
from pystencils.field import Field

with CodeGeneration() as ctx:
    streaming_pattern = 'aa'
    target = 'cpu'
    stencil = get_stencil('D3Q19')
    dim = len(stencil[0])
    values_per_cell = len(stencil)
    collision_rule = create_lb_collision_rule(method='srt', stencil=stencil)
    pdf_field = Field.create_generic('pdfs', dim, index_shape=(values_per_cell,), layout='fzyx')
    optimization = {
        'symbolic_field': pdf_field,
        'target': target
    }

    #   Generate PackInfo specifically for streaming pattern
    generate_lb_pack_info(ctx, 'AccessorBasedPackInfo', stencil, pdf_field,
                          streaming_pattern=streaming_pattern, target=target, namespace='pystencils')

    #   Generate reference using the alternating pull/push approach
    update_rule_odd = create_lb_update_rule(collision_rule=collision_rule, optimization=optimization,
                                            streaming_pattern=streaming_pattern, timestep=Timestep.ODD)
    generate_pack_info_from_kernel(ctx, 'FromKernelPackInfoPull', update_rule_odd, kind='pull', target=target)
    generate_pack_info_from_kernel(ctx, 'FromKernelPackInfoPush', update_rule_odd, kind='push', target=target)
