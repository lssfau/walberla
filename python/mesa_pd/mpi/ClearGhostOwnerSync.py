# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class ClearGhostOwnerSync:
    def __init__(self):
        self.context = {'properties': [], 'interface': []}

        self.context['interface'].append(
            create_access("flags", 'walberla::mesa_pd::data::particle_flags::FlagT', access='g'))
        self.context['interface'].append(create_access('ghostOwners', 'std::vector<int>', access='r'))
        self.context['interface'].append(create_access('neighborState', 'std::unordered_set<walberla::mpi::MPIRank>',
                                                       access='r'))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'mpi/ClearGhostOwnerSync.templ.h', ctx)
