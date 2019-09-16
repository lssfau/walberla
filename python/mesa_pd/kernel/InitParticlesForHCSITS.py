# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from ..Container import Container
from mesa_pd.utility import generateFile

class InitParticlesForHCSITS(Container):
   def __init__(self):
      super().__init__()
      self.addProperty("globalAcceleration", "walberla::mesa_pd::Vec3", defValue = "0")

      self.accessor = Accessor()
      self.accessor.require("uid",             "walberla::id_t",                                    access="g")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3",                           access="gr")
      self.accessor.require("angularVelocity", "walberla::mesa_pd::Vec3",                           access="gr")
      self.accessor.require("invMass",         "walberla::real_t",                                  access="g" )
      self.accessor.require("inertia",         "walberla::mesa_pd::Mat3",                           access="g" )
      self.accessor.require("invInertia",      "walberla::mesa_pd::Mat3",                           access="g" )
      self.accessor.require("dv",              "walberla::mesa_pd::Vec3",                           access="gr" )
      self.accessor.require("dw",              "walberla::mesa_pd::Vec3",                           access="gr" )
      self.accessor.require("torque",          "walberla::mesa_pd::Vec3",                           access="r" )
      self.accessor.require("force",           "walberla::mesa_pd::Vec3",                           access="r" )

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["properties"]      = self.properties
      context["interface"]       = self.accessor.properties
      generateFile(path, 'kernel/InitParticlesForHCSITS.templ.h', context)
