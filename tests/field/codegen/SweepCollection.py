import sympy as sp

import pystencils as ps
from pystencils import Assignment
from pystencils_walberla import CodeGeneration, function_generator, generate_sweep_collection


with CodeGeneration() as ctx:
    field_type = "float64" if ctx.double_accuracy else "float32"

    a = sp.Symbol('a')
    f1, f2, f3 = ps.fields(f"f1, f2, f3: {field_type}[3D]", layout='fzyx')
    up1 = Assignment(f2.center, 2 * a * f1.center)
    up2 = Assignment(f3.center, 2 * a * f2.center)

    fct1 = function_generator(ctx, 'fct1', up1)
    fct2 = function_generator(ctx, 'fct2', up2)

    generate_sweep_collection(ctx, "SweepCollection", [fct1, fct2])
