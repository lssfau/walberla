import pystencils as ps
from lbmpy.updatekernels import create_stream_pull_only_kernel
from lbmpy.stencils import get_stencil
from pystencils_walberla import CodeGeneration, generate_sweep

with CodeGeneration() as ctx:
    f_size = 19
    dtype = 'float64' if ctx.double_accuracy else 'float32'

    # Copy sweep
    src, dst = ps.fields("src({f_size}), dst({f_size}) : {dtype}[3D]".format(dtype=dtype, f_size=f_size),
                         layout='fzyx')
    copy_only = [ps.Assignment(dst(i), src(i)) for i in range(f_size)]
    generate_sweep(ctx, 'MicroBenchmarkCopyKernel', copy_only,
                   target='gpu', gpu_indexing_params={'block_size': (128, 1, 1)})

    # Stream-only sweep
    stencil = get_stencil("D3Q19")
    stream_only = create_stream_pull_only_kernel(stencil, src_field_name='src', dst_field_name='dst',
                                                 generic_field_type=dtype, generic_layout='fzyx')
    generate_sweep(ctx, 'MicroBenchmarkStreamKernel', stream_only,
                   target='gpu', gpu_indexing_params={'block_size': (128, 1, 1)})
