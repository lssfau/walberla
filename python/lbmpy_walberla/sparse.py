import numpy as np

from lbmpy.sparse import (
    create_lb_update_rule_sparse, create_macroscopic_value_getter_sparse,
    create_macroscopic_value_setter_sparse, create_symbolic_list)
from pystencils import TypedSymbol, create_kernel


class ListLbGenerator:
    def __init__(self, collision_rule, pdf_type=np.float64, index_type=np.uint32, layout='SoA'):
        self.collision_rule = collision_rule

        num_cells = TypedSymbol('num_cells', np.uint64)
        method = self.collision_rule.method
        q = len(method.stencil)

        self.num_cells = num_cells
        self.src = create_symbolic_list('src', num_cells, q, pdf_type, layout)
        self.dst = create_symbolic_list('dst', num_cells, q, pdf_type, layout)
        self.idx = create_symbolic_list('idx', num_cells, q, index_type, layout)

    def kernel(self, kernel_type='stream_pull_collide'):
        ac = create_lb_update_rule_sparse(self.collision_rule, self.src, self.dst, self.idx, kernel_type)
        return create_kernel(ac)

    def getter_ast(self, density=True, velocity=True):
        assert density or velocity

        method = self.collision_rule.method
        output_descriptor = {}
        if density:
            output_descriptor['density'] = create_symbolic_list('rho', self.src.spatial_shape[0], 1, self.src.dtype)
        if velocity:
            output_descriptor['velocity'] = create_symbolic_list('vel', self.src.spatial_shape[0],
                                                                 method.dim, self.src.dtype)

        ac = create_macroscopic_value_getter_sparse(method, self.src, output_descriptor)
        return create_kernel(ac)

    def setter_ast(self, density=True, velocity=True):
        method = self.collision_rule.method
        if density is True:
            density = create_symbolic_list('rho', self.src.spatial_shape[0], 1, self.src.dtype).center

        if velocity is True:
            velocity = create_symbolic_list('vel', self.src.spatial_shape[0], method.dim, self.src.dtype).center_vector

        ac = create_macroscopic_value_setter_sparse(method, self.src, density, velocity)
        return create_kernel(ac)
