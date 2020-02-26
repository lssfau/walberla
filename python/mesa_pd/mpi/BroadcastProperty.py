# -*- coding: utf-8 -*-

from ..utility import generate_file


class BroadcastProperty:
    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'mpi/BroadcastProperty.templ.h', ctx)
