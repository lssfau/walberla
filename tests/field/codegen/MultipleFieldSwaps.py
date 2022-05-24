import pystencils as ps
from pystencils_walberla import CodeGeneration, generate_sweep

with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"
    src_1, dst_1 = ps.fields(f"src1, src1_tmp: {field_type}[2D]", layout='fzyx')
    src_2, dst_2 = ps.fields(f"src2, src2_tmp: {field_type}[2D]", layout='fzyx')
    src_3, dst_3 = ps.fields(f"src3, src3_tmp: {field_type}[2D]", layout='fzyx')

    assignments = ps.AssignmentCollection([ps.Assignment(dst_1.center, 2 * src_1.center),
                                           ps.Assignment(dst_2.center, 2 * src_2.center),
                                           ps.Assignment(dst_3.center, 2 * src_3.center)])

    generate_sweep(ctx, 'MultipleFieldSwaps', assignments,
                   field_swaps=[(src_1, dst_1), (src_2, dst_2), (src_3, dst_3)])
