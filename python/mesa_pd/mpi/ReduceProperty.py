# -*- coding: utf-8 -*-

from ..utility import generate_file


class ReduceProperty:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'mpi/ReduceProperty.templ.h', ctx)
