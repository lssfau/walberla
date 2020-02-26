# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class InsertParticleIntoLinkedCells:
    def __init__(self):
        self.context = {'interface': []}
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(
            create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g"))
        self.context['interface'].append(create_access("nextParticle", "size_t", access="gs"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/InsertParticleIntoLinkedCells.templ.h', ctx)
