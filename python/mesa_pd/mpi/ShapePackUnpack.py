# -*- coding: utf-8 -*-

from ..utility import generate_file


class ShapePackUnpack:
    def __init__(self, particle_storage):
        self.ps = particle_storage

    def generate(self, module):
        ctx = {'module': module, 'particle': self.ps.get_context()}
        generate_file(module['module_path'], 'mpi/ShapePackUnpack.templ.h', ctx)
