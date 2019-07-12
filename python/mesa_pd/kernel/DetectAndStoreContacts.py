# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class DetectAndStoreContacts:

   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("uid",             "walberla::id_t",                                    access="g")
      self.accessor.require("flags",           "walberla::mesa_pd::data::particle_flags::FlagT",    access="g")
      self.accessor.require("position",        "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("rotation",        "walberla::mesa_pd::Rot3",                           access="g")
      self.accessor.require("shape",           "BaseShape*",                                        access="g")



   def getRequirements(self):
      return self.accessor


   def generate(self, path):
      context = dict()
      context["interface"]        = self.accessor.properties
      generateFile(path, 'kernel/DetectAndStoreContacts.templ.h', context)
