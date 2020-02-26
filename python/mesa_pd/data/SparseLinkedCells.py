# -*- coding: utf-8 -*-

import numpy as np
from ..utility import generate_file


class SparseLinkedCells:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'data/SparseLinkedCells.templ.h', ctx)
