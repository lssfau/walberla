import pystencils as ps
from lbmpy.updatekernels import create_stream_only_kernel
from lbmpy import LBStencil, Stencil
from pystencils_walberla import CodeGeneration, generate_sweep

with CodeGeneration() as ctx:
    f_size = 19
    dtype = 'float64' if ctx.double_accuracy else 'float32'

    # Copy sweep
    src, dst = ps.fields(f"src({f_size}), dst({f_size}) : {dtype}[3D]", layout='fzyx')

    copy_only = [ps.Assignment(dst(i), src(i)) for i in range(f_size)]
    generate_sweep(ctx, 'MicroBenchmarkCopyKernel', copy_only,
                   target=ps.Target.GPU, gpu_indexing_params={'block_size': (128, 1, 1)})

    # Stream-only sweep
    stencil = LBStencil(Stencil.D3Q19)
    stream_only = create_stream_only_kernel(stencil, src_field=src, dst_field=dst)
    generate_sweep(ctx, 'MicroBenchmarkStreamKernel', stream_only,
                   target=ps.Target.GPU, gpu_indexing_params={'block_size': (128, 1, 1)})
