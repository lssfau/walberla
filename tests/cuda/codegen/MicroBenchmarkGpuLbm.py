import pystencils as ps
from pystencils_walberla.sweep import Sweep
from lbmpy.updatekernels import create_stream_pull_only_kernel
from lbmpy.stencils import get_stencil

dtype = 'float64'
f_size = 19


def copy_only():
    src, dst = ps.fields("src({f_size}), dst({f_size}) : {dtype}[3D]".format(dtype=dtype, f_size=f_size),
                         layout='fzyx')
    return [ps.Assignment(dst(i), src(i)) for i in range(f_size)]


def stream_only():
    stencil = get_stencil("D3Q19")
    return create_stream_pull_only_kernel(stencil, src_field_name='src',
                                          dst_field_name='dst',
                                          generic_field_type=dtype,
                                          generic_layout='fzyx')


opt = {'gpu_indexing_params': {'block_size': (128, 1, 1)}, 'data_type': dtype}

Sweep.generate_from_equations('MicroBenchmarkCopyKernel', copy_only, target='gpu', optimization=opt)
Sweep.generate_from_equations('MicroBenchmarkStreamKernel', stream_only, target='gpu', optimization=opt)
