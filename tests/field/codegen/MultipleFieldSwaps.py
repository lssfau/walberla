import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep

with CodeGeneration() as ctx:
    src_1, dst_1 = ps.fields("src1, src1_tmp: [2D]", layout='fzyx')
    src_2, dst_2 = ps.fields("src2, src2_tmp: [2D]", layout='fzyx')
    src_3, dst_3 = ps.fields("src3, src3_tmp: [2D]", layout='fzyx')

    assignments = ps.AssignmentCollection([ps.Assignment(dst_1.center, 2 * src_1.center),
                                           ps.Assignment(dst_2.center, 2 * src_2.center),
                                           ps.Assignment(dst_3.center, 2 * src_3.center)])

    generate_sweep(ctx, 'MultipleFieldSwaps', assignments,
                   field_swaps=[(src_1, dst_1), (src_2, dst_2), (src_3, dst_3)])
