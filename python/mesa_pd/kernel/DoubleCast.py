# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class DoubleCast:
    def __init__(self, particle_storage):
        self.ps = particle_storage
        self.context = {'interface': [create_access("shape", "BaseShape*", access="g")]}

    def generate(self, module):
        context = {'module': module, 'particle': self.ps.get_context(), **self.context}
        generate_file(module['module_path'], 'kernel/DoubleCast.templ.h', context)
