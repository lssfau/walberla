# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class ClearNextNeighborSync:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("flags",        "walberla::mesa_pd::data::particle_flags::FlagT", access="g")
      self.accessor.require("ghostOwners",  "std::vector<int>",                           access="r")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["interface"] = self.accessor.properties
      generateFile(path, 'mpi/ClearNextNeighborSync.templ.h', context)
