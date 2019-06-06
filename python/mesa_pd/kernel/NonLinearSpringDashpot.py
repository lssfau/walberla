# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class NonLinearSpringDashpot:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("uid",             "walberla::id_t",                                access="g")
      self.accessor.require("position",        "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("force",           "walberla::mesa_pd::Vec3",                           access="r" )
      self.accessor.require("angularVelocity", "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("torque",          "walberla::mesa_pd::Vec3",                           access="r" )
      self.accessor.require("type",            "uint_t",                                        access="g")
      self.accessor.require("contactHistory",  "std::map<walberla::id_t, walberla::mesa_pd::Vec3>", access="gs")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["parameters"]       = ["lnCORsqr", "meff", "stiffnessT", "dampingT", "frictionCoefficientStatic", "frictionCoefficientDynamic"]
      context["interface"]        = self.accessor.properties

      generateFile(path, 'kernel/NonLinearSpringDashpot.templ.h', context)
