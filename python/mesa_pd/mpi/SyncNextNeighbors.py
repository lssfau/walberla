# -*- coding: utf-8 -*-

from ..utility import generate_file


class SyncNextNeighbors:
    def __init__(self, particle_storage):
        particle_storage.add_property("position", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
        particle_storage.add_property("interactionRadius", "walberla::real_t", defValue="real_t(0)",
                                      syncMode="ON_GHOST_CREATION")
        particle_storage.add_property("flags", "walberla::mesa_pd::data::particle_flags::FlagT", defValue="",
                                      syncMode="ON_GHOST_CREATION")
        particle_storage.add_property("owner", "int", defValue="-1", syncMode="ON_GHOST_CREATION")
        particle_storage.add_property("ghostOwners", "std::unordered_set<walberla::mpi::MPIRank>", defValue="",
                                      syncMode="NEVER")

    def generate(self, module):
        ctx = {'module': module}
        generate_file(module['module_path'], 'mpi/SyncNextNeighbors.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/SyncNextNeighbors.templ.cpp', ctx)
