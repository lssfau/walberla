# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from ..Container import Container
from mesa_pd.utility import generateFile

class HCSITSRelaxationStep(Container):
   def __init__(self):
      super().__init__()
      self.addProperty("maxSubIterations", "size_t", defValue = "20")
      self.addProperty("relaxationModel", "RelaxationModel", defValue = "InelasticFrictionlessContact")
      self.addProperty("deltaMax", "real_t", defValue = "0")
      self.addProperty("cor", "real_t", defValue = "real_t(0.2)")

      self.accessor = Accessor()
      self.accessor.require("uid",             "walberla::id_t",                                    access="g")
      self.accessor.require("position",        "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("angularVelocity", "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("invMass",         "walberla::real_t",                                  access="g" )
      self.accessor.require("invInertia",      "walberla::mesa_pd::Mat3",                           access="g" )
      self.accessor.require("dv",              "walberla::mesa_pd::Vec3",                           access="gr" )
      self.accessor.require("dw",              "walberla::mesa_pd::Vec3",                           access="gr" )

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["properties"]      = self.properties
      context["interface"]        = self.accessor.properties
      generateFile(path, 'kernel/HCSITSRelaxationStep.templ.h', context)
