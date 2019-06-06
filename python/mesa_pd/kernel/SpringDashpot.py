# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile, checkInterface

class SpringDashpot:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("position",        "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("force",           "walberla::mesa_pd::Vec3",                           access="r" )
      self.accessor.require("angularVelocity", "walberla::mesa_pd::Vec3",                           access="g")
      self.accessor.require("torque",          "walberla::mesa_pd::Vec3",                           access="r" )
      self.accessor.require("type",            "uint_t",                                            access="g")
      self.accessor.require("contactHistory",  "std::map<walberla::id_t, walberla::mesa_pd::Vec3>", access="gs")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["parameters"]       = ["stiffness", "dampingN", "dampingT", "friction"]
      context["interface"]        = self.accessor.properties

      generateFile(path, 'kernel/SpringDashpot.templ.h', context)

      context["InterfaceTestName"] = "SpringDashpotInterfaceCheck"
      context["KernelInclude"] = "kernel/SpringDashpot.h"
      context["ExplicitInstantiation"] = "template void kernel::SpringDashpot::operator()(const size_t p_idx1, const size_t p_idx2, Accessor& ac, const Vec3& contactPoint, const Vec3& contactNormal, const real_t& penetrationDepth) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/SpringDashpotInterfaceCheck.cpp')
