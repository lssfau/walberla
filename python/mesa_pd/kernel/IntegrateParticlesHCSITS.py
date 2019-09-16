# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from ..Container import Container
from mesa_pd.utility import generateFile

class IntegrateParticlesHCSITS(Container):
   def __init__(self):
      super().__init__()
      self.addProperty("speedLimiterActive", "bool", defValue = "false")
      self.addProperty("speedLimitFactor",    "real_t", defValue ="real_t(1.0)")

      self.accessor = Accessor()
      self.accessor.require("uid",              "walberla::id_t",                                   access="g")
      self.accessor.require("position",         "walberla::mesa_pd::Vec3",                          access="gr")
      self.accessor.require("rotation",         "walberla::mesa_pd::Rot3",                          access="gr")
      self.accessor.require("linearVelocity",   "walberla::mesa_pd::Vec3",                          access="gr")
      self.accessor.require("angularVelocity",  "walberla::mesa_pd::Vec3",                          access="gr")
      self.accessor.require("dv",               "walberla::mesa_pd::Vec3",                          access="g" )
      self.accessor.require("dw",               "walberla::mesa_pd::Vec3",                          access="g" )
      self.accessor.require("flags",            "walberla::mesa_pd::data::particle_flags::FlagT",   access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["properties"]      = self.properties
      context["interface"]       = self.accessor.properties
      generateFile(path, 'kernel/IntegrateParticlesHCSITS.templ.h', context)
