# -*- coding: utf-8 -*-

from ..utility import generate_file


class SyncNextNeighborsNoGhosts:
    def __init__(self, particle_storage):
        particle_storage.add_property("position", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
        particle_storage.add_property("owner", "int", defValue="-1", syncMode="ON_GHOST_CREATION")

    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'mpi/SyncNextNeighborsNoGhosts.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/SyncNextNeighborsNoGhosts.templ.cpp', ctx)
