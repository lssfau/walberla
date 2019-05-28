# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile, checkInterface

class HeatConduction:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("temperature",     "walberla::real_t", access="g")
      self.accessor.require("heatFlux",        "walberla::real_t", access="gsr")
      self.accessor.require("type",            "uint_t",           access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["parameters"]       = ["conductance"]
      context["interface"]        = self.accessor.properties

      generateFile(path, 'kernel/HeatConduction.templ.h', context)

      context["InterfaceTestName"] = "HeatConductionInterfaceCheck"
      context["KernelInclude"] = "kernel/HeatConduction.h"
      context["ExplicitInstantiation"] = "template void kernel::HeatConduction::operator()(const size_t p_idx1, const size_t p_idx2, Accessor& ac) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/HeatConductionInterfaceCheck.cpp')
