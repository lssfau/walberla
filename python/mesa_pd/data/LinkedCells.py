# -*- coding: utf-8 -*-

import numpy as np
from ..utility import generate_file


class LinkedCells:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'data/LinkedCells.templ.h', ctx)
