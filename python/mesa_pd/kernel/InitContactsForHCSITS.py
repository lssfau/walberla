# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from ..Container import Container
from mesa_pd.utility import generateFile

class InitContactsForHCSITS(Container):
   def __init__(self):
      super().__init__()
      self.addProperty("erp", "real_t", defValue = "real_t(0.8)")
      self.addProperty("maximumPenetration", "real_t", defValue ="0")

      self.paccessor = Accessor()
      self.paccessor.require("uid",             "walberla::id_t",                                    access="g")
      self.paccessor.require("position",        "walberla::mesa_pd::Vec3",                           access="g")
      self.paccessor.require("invInertia",      "walberla::mesa_pd::Mat3",                           access="g" )
      self.paccessor.require("invMass",         "walberla::real_t",                                  access="g" )


   def getRequirements(self):
      return self.paccessor

   def generate(self, path):
      context = dict()
      context["properties"]      = self.properties
      context["material_parameters"] = ["friction"]
      context["interface"]        = self.paccessor.properties
      generateFile(path, 'kernel/InitContactsForHCSITS.templ.h', context)
