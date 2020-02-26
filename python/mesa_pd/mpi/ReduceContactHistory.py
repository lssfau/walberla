# -*- coding: utf-8 -*-

from ..utility import generate_file


class ReduceContactHistory:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'mpi/ReduceContactHistory.templ.h', ctx)
