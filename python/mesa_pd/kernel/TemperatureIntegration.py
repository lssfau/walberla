# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class TemperatureIntegration:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("temperature",     "walberla::real_t", access="gs")
      self.accessor.require("heatFlux",        "walberla::real_t", access="gs")
      self.accessor.require("type",            "uint_t",           access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["parameters"]       = ["invHeatCapacity"]
      context["interface"] = self.accessor.properties
      generateFile(path, 'kernel/TemperatureIntegration.templ.h', context)

      context["InterfaceTestName"] = "TemperatureIntegrationInterfaceCheck"
      context["KernelInclude"] = "kernel/TemperatureIntegration.h"
      context["ExplicitInstantiation"] = "template void kernel::TemperatureIntegration::operator()(const size_t p_idx1, Accessor& ac) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/TemperatureIntegrationInterfaceCheck.cpp')
