import pystencils as ps

from lbmpy.creationfunctions import create_lb_update_rule, LBMConfig, LBMOptimisation
from lbmpy.advanced_streaming import Timestep
from lbmpy.enums import Method, Stencil
from lbmpy.stencils import LBStencil
from pystencils_walberla import CodeGeneration, generate_pack_info_from_kernel
from lbmpy_walberla.packinfo import generate_lb_pack_info
from pystencils.field import Field

with CodeGeneration() as ctx:
    data_type = "float64" if ctx.double_accuracy else "float32"
    streaming_pattern = 'aa'
    target = ps.Target.CPU
    stencil = LBStencil(Stencil.D3Q19)

    pdf_field = Field.create_generic('pdfs', stencil.D, index_shape=(stencil.Q,), layout='fzyx', dtype=data_type)

    lbm_config = LBMConfig(method=Method.SRT, stencil=stencil, streaming_pattern=streaming_pattern,
                           timestep=Timestep.ODD)
    lbm_opt = LBMOptimisation(symbolic_field=pdf_field)

    #   Generate PackInfo specifically for streaming pattern
    generate_lb_pack_info(ctx, 'AccessorBasedPackInfo', stencil, pdf_field,
                          streaming_pattern=streaming_pattern, target=target, namespace='pystencils')

    #   Generate reference using the alternating pull/push approach
    update_rule_odd = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    generate_pack_info_from_kernel(ctx, 'FromKernelPackInfoPull', update_rule_odd, kind='pull', target=target)
    generate_pack_info_from_kernel(ctx, 'FromKernelPackInfoPush', update_rule_odd, kind='push', target=target)
