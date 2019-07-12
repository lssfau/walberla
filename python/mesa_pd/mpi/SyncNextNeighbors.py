# -*- coding: utf-8 -*-

from ..utility import generateFile

class SyncNextNeighbors:
   def __init__(self, p):
      p.addProperty("position",          "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
      p.addProperty("interactionRadius", "walberla::real_t",        defValue="real_t(0)", syncMode="ONCE")
      p.addProperty("flags",             "walberla::mesa_pd::data::particle_flags::FlagT", defValue="", syncMode="ONCE")
      p.addProperty("owner",             "int",                     defValue="-1",        syncMode="ONCE")
      p.addProperty("ghostOwners",       "std::unordered_set<walberla::mpi::MPIRank>",    defValue="",          syncMode="NEVER")

   def generate(self, path):
      generateFile(path, 'mpi/SyncNextNeighbors.templ.h')
      generateFile(path, 'mpi/SyncNextNeighbors.templ.cpp')
