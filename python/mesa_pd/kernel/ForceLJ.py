# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class ForceLJ:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("position",        "walberla::mesa_pd::Vec3", access="g")
      self.accessor.require("force",           "walberla::mesa_pd::Vec3", access="r" )
      self.accessor.require("type",            "uint_t",              access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["parameters"]       = ["epsilon", "sigma"]
      context["interface"]        = self.accessor.properties

      generateFile(path, 'kernel/ForceLJ.templ.h', context)

      context["InterfaceTestName"] = "ForceLJInterfaceCheck"
      context["KernelInclude"] = "kernel/ForceLJ.h"
      context["ExplicitInstantiation"] = "template void kernel::ForceLJ::operator()(const size_t p_idx1, const size_t p_idx2, Accessor& ac) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/ForceLJInterfaceCheck.cpp')
