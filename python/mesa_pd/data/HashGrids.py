# -*- coding: utf-8 -*-

from ..utility import generate_file


class HashGrids:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'data/HashGrids.templ.h', ctx)
