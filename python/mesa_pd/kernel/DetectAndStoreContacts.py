# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class DetectAndStoreContacts:

    def __init__(self):
        self.context = {
            'interface': [
                create_access("uid", "walberla::id_t", access="g"),
                create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g"),
                create_access("position", "walberla::mesa_pd::Vec3", access="g"),
                create_access("rotation", "walberla::mesa_pd::Rot3", access="g"),
                create_access("shape", "BaseShape*", access="g")
            ]
        }

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/DetectAndStoreContacts.templ.h', ctx)
